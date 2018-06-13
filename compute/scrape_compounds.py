import os
import time
import re
import urllib.request

import pubchempy as pcp

from bs4 import BeautifulSoup

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "June 2018"

# Notes:
# record_type should be 3d, if possible (2d by default)
# Could search by superstructure (default options should be good, except maybe ChainsMatchRings - defaults true)
# Get SMILES string of cyclohexene, use that as superstructure
# Could search by similarity (play around with Threshold - default 90)

cids = {"products": [8079],
        "dienes": [7845],
        "dienophiles": [6325]
        }

#TODO: Figure out it PUG will give 3D coordinates by default


class PUGScraper:
    """
    Uses the PubChem Power User Gateway (PUG) as well as the PUG REST API to
    query and identify molecules similar to a set of base molecules. In
    addition to scraping chemical ids, this class can download structure files
    and process them with pybel.

    Thank you to http://python.zirael.org/e-pubchem1.html for some inspiration
    in devising this class.
    """

    def __init__(self, cids, base_dir, sub_dirs=None, use_rest=False):
        """
        PUGScraper
        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest
        :param base_dir: A string representing an absolute path where downloaded
        information (ids, SDF files) should go.
        :param sub_dirs: By default, this is None. Otherwise, it should be a
        dict {"category": "sub_dir"}, where each category in cids is
        represented, and where each sub_dir is a string to a subdirectory within
        base_dir (these will be created if they don't already exist).
        :param use_rest: If True, prefer to use the PUG-REST API through
        pubchempy. By default, this is False, and the full PUG will be used,
        to avoid timeouts. It is preferred that
        """

        self.cids = cids
        self.base_dir = base_dir
        self.sub_dirs = sub_dirs
        self.use_rest = use_rest

        if not self.use_rest:
            # We only need to care about URLs if we're using PUG
            # For the REST API, pubchempy will handle it
            self.base_url = "http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
            self.cid_url = '''http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=xml&rettype=uilist&WebEnvRq=1&db=pccompound&query_key=%(key)s&WebEnv=%(env)s'''

            # Since we only care about a small subset of possible job types,
            # simple string formatting will suffice.
            # TODO: consider using a template engine, or dedicated XML library
            self.template_query = """<?xml version="1.0"?>
<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
<PCT-Data>
  <PCT-Data_input>
    <PCT-InputData>
      <PCT-InputData_query>
        <PCT-Query>
          <PCT-Query_type>
            <PCT-QueryType>
              <PCT-QueryType_css>
                <PCT-QueryCompoundCS>
                  <PCT-QueryCompoundCS_query>
                    <PCT-QueryCompoundCS_query_data>%(id)s</PCT-QueryCompoundCS_query_data>
                  </PCT-QueryCompoundCS_query>
                  <PCT-QueryCompoundCS_type>
                    <PCT-QueryCompoundCS_type_%(query_type)s>
                      %(type_specific)s
                    </PCT-QueryCompoundCS_type_%(query_type)s>
                  </PCT-QueryCompoundCS_type>
                  <PCT-QueryCompoundCS_results>%(max_records)s</PCT-QueryCompoundCS_results>
                </PCT-QueryCompoundCS>
              </PCT-QueryType_css>
            </PCT-QueryType>
          </PCT-Query_type>
        </PCT-Query>
      </PCT-InputData_query>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data> 
"""

            self.template_download = """<?xml version="1.0"?>
<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
<PCT-Data>
  <PCT-Data_input>
    <PCT-InputData>
      <PCT-InputData_download>
        <PCT-Download>
          <PCT-Download_uids>
            <PCT-QueryUids>
              <PCT-QueryUids_ids>
                <PCT-ID-List>
                  <PCT-ID-List_db>pccompound</PCT-ID-List_db>
                  <PCT-ID-List_uids>
                    %(uid_strings)s
                  </PCT-ID-List_uids>
                </PCT-ID-List>
              </PCT-QueryUids_ids>
            </PCT-QueryUids>
          </PCT-Download_uids>
          <PCT-Download_format value="%(format)s"/>
          <PCT-Download_compression value="gzip"/>
        </PCT-Download>
      </PCT-InputData_download>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data>"""

            self.template_check = """
<PCT-Data>
  <PCT-Data_input>
    <PCT-InputData>
      <PCT-InputData_request>
        <PCT-Request>
          <PCT-Request_reqid>%(request_id)s</PCT-Request_reqid>
          <PCT-Request_type value="status"/>
        </PCT-Request>
      </PCT-InputData_request>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data>            
"""

    def create_dirs(self):
        """
        Ensures that all paths and subpaths required by PUGScraper exist.
        :return:
        """
        try:
            if not os.path.exists(self.base_dir):
                os.makedirs(directory)
            if self.sub_dirs is not None:
                for cat in self.sub_dirs.keys():
                    sub_dir = os.path.join(self.base_dir, self.sub_dirs[cat])
                    if not os.path.exists(sub_dir):
                        os.makedirs(sub_dir)
        except:
            raise ValueError("Given directory paths are invalid; check that"
                             "strings have been input correctly, and that you"
                             "have permission to use the given directories.")

    def check_queries(self, queries):
        if queries == 5:
            time.sleep(1)
            return 0
        else:
            return queries + 1

    def check_and_process(self, output, path, request_type):
        """
        Parses output XML. If the output returns a file, downloads that file
        at the specified path; if the output returns a request id, continue to
        check at that id until the request is complete.

        :param output: An XML file generated by download_files_pug
        :param path: A path string for downloads, generated by
        download_files_pug
        :return: If request_type is "query", return list of CIDs; if
        request_type is "download", return nothing.
        """

        wait_time = 2
        num_checks = 450

        for i in range(num_checks):
            parsed = BeautifulSoup(output)

            if request_type == "download":
                finished = parsed.find("PCT-Download-URL_url")

                if finished is not None:
                    ftp_url = finished.text
                    opener = urllib.request.URLopener()
                    opener.retrieve(ftp_url, path)
                    return

            elif request_type == "query":
                query = parsed.find("PCT-Entrez_query-key")
                env = parsed.find("PCT-Entrez_webenv")

                if query is not None and env is not None:
                    query = query.text
                    env = env.text

                    url = self.cid_url % {"query": query, "env": env}

                    results = urllib.request.urlopen(url)
                    cid_output = results.read()
                    results.close()

                    cids_parsed = BeautifulSoup(cid_output)
                    return [int(id) for id in cids_parsed.findAll("Id")]


            time.sleep(wait_time)

            request_id = parsed.find("PCT-Request_reqid").text
            request = self.template_check % {"request_id": request_id}
            page = urllib.request.urlopen(self.base_url, data=request)
            output = page.read()
            page.close()

    def download_files(self, cids, pngs=True, download_parents=False):
        """
        Generalized function for downloading files (both SDF and PNG, for quick
        reference of structure and for full coordinate and bonding information),
        which calls either download_files_rest or download_files_pug, depending
        on if REST is being used.

        :param cids: A dict {"category": {id:[ids]}}, where each category is a
        molecular type of interest.
        :param pngs: If True, PNG files will be downloaded alongside SDF files.
        :param download_parents: If True, then files will be downloaded for
        parent molecules, in addition to the molecules returned from their
        queries.
        :return:
        """

        self.create_dirs()

        if self.use_rest:
            self.download_files_rest(cids, pngs=pngs,
                                     download_parents=download_parents)
        else:
            self.download_files_pug(cids, pngs=pngs,
                                    download_parents=download_parents)

    def download_files_rest(self, cids, pngs=True, download_parents=False):
        """
        Downloads SDF files (and, potentially, PNG files) from PUG-REST.

        The PUB-REST API only allows a maximum of 5 requests per second and 400
        total requests per minute. In order to ensure that this function never
        exceeds this limit, after every five queries, the program will halt for
        one second.
        TODO: Find a more elegant way to handle this

        NOTE: Each compound will receive a file cid.sdf if it is a "parent"
        molecule (it was used to start a query), where id is the CID of the
        molecule, or id_parent if it was a hit from a query, where id is the CID
        of the molecule, and parent is the CID from the parent molecule. This
        naming is a way to keep track of where each file came from.
        TODO: Is this a good naming convention? How should I go about naming?

        :param cids: A dict {"category": {id:[ids]}}, where each category is a
        molecular type of interest.
        :param pngs: If True, PNG files will be downloaded alongside SDF files.
        :return:
        """

        queries_run = 0

        formats = ["SDF"]
        if pngs:
            formats.append("PNG")

        for cat in cids.keys():
            if self.sub_dirs is not None:
                cat_path = os.path.join(self.base_dir, self.sub_dirs[cat])
            else:
                cat_path = os.path.join(self.base_dir)

            for parent in cids[cat].keys():

                if download_parents:
                    for format in formats:
                        filename = str(parent) + "." + format.lower()
                        filepath = os.path.join(cat_path, filename)

                        queries_run = self.check_queries(queries_run)
                        pcp.download(format, filepath, parent)

                for cid in cids[cat][parent]:
                    for format in formats:
                        filename = str(cid) + "_" + str(parent) + "." + format.lower()
                        filepath = os.path.join(cat_path, filename)

                        queries_run = self.check_queries(queries_run)
                        pcp.download(format, filepath, cid)

    def download_files_pug(self, cids, pngs=True, download_parents=False):
        """
        Downloads SDF files (and, potentially, PNG files) from PUG.

        NOTE: PUG will zip the results of each job together. This means that we
        cannot add our own naming scheme to each file, as in
        download_files_rest. However, we will still split the job into
        categories. This may be somewhat inefficient, but helps with data
        processing.

        :param cids: A dict {"category": {id:[ids]}}, where each category is a
        molecular type of interest.
        :param pngs: If True, PNG files will be downloaded alongside SDF files.
        :return:
        """

        formats = ["sdf"]
        if pngs:
            formats.append("png")

        uid_template = "<PCT-ID-List_uids_E>%(cid)s</PCT-ID-List_uids_E>"

        for cat in cids.keys():
            download_ids = []

            if self.sub_dirs is not None:
                cat_path = os.path.join(self.base_dir, self.sub_dirs[cat])
            else:
                cat_path = os.path.join(self.base_dir)

            for parent in cids[cat].keys():

                if download_parents:
                    download_ids.append(parent)

                for cid in cids[cat][parent]:
                    download_ids.append(cid)

            uid_templates = [uid_template % {"cid": str(cid)}
                             for cid in download_ids]
            uid_strings = "\n".join(uid_templates)

            for format in formats:
                xml = self.template_download % {"format":format,
                                                "uid_strings": uid_strings}
                page = urllib.request.urlopen(self.base_url, data=xml)
                output = page.read()
                page.close()
                filename = str(cat) + "_." + format + ".gz"
                filepath = os.path.join(cat_path, filename)
                self.check_and_process(output, filepath)

    def store_cids(self, cids):
        """
        Store cids in a text (.txt) file.

        By default, this stores not only the category and the cids, but also the
        CID that caused each other CID to be matched.

        :param cids: A dict {"category": {id:[ids]}}, where each category is a
        molecular type of interest.
        :return:
        """

        self.create_dirs()

        for cat in cids.keys():
            if self.sub_dirs is not None:
                filepath = os.path.join(self.base_dir, self.sub_dirs[cat],
                                        "cids.txt")
            else:
                filepath = os.path.join(self.base_dir, "cids.txt")

            with open(filepath, 'a') as cidfile:
                cidfile.write("Category: " + cat + "\n")
                for parent in cids[cat].keys():
                    cidfile.write("\t" + "Parent: " + str(parent) + "\n")
                    for cid in cids[cat][parent]:
                        cidfile.write("\t\t" + str(cid) + "\n")

    def scrape_similar(self, cids, threshold=90, max_records=10000):
        """
        Generalized function for similarity searches (which look for molecules
        similar to the query) that either calls scrape_similar_rest
        scrape_similar_pug, depending on if REST is being used. Kwargs are
        those used by PUG and PUG-REST for similarity queries.

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :param threshold: int representing the minimum Tanimoto score required
        for a hit.
        :param max_records: int representing the maximum number of hits allowed.
        NOTE: for the REST API, this should be kept relatively low to avoid
        timeout.
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        if self.use_rest:
            self.scrape_similar_rest(cids, threshold=threshold,
                                     max_records=max_records)
        else:
            self.scrape_similar_pug(cids, threshold=threshold,
                                    max_records=max_records)

    def scrape_similar_rest(self, cids, threshold=90, max_records=10000):
        """
        Searches by similarity (2D) using PUG-REST.

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :param threshold: int representing the minimum Tanimoto score required
        for a hit.
        :param max_records: int representing the maximum number of hits allowed.
        NOTE: for the REST API, this should be kept relatively low to avoid
        timeout.
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        output = {}
        queries_run = 0

        for cat in cids.keys():
            output[cat] = {}
            for cid in cids[cat]:
                queries_run = self.check_queries(queries_run)
                result = pcp.get_cids(cid, searchtype="similar",
                                      threshold=threshold,
                                      max_records=max_records, record_type="3d")
                output[cat][cid] = result

        return output

    def scrape_similar_pug(self, cids, threshold=90, max_records=10000):
        """
        Searches by similarity (2D) using PUG.

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :param threshold: int representing the minimum Tanimoto score required
        for a hit.
        :param max_records: int representing the maximum number of hits allowed.
        NOTE: for the REST API, this should be kept relatively low to avoid
        timeout.
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        cid_results = {}

        threshold_tag = """"<PCT-CSSimilarity_threshold>
%(threshold)s
</PCT-CSSimilarity_threshold>""" % {"threshold": str(threshold)}

        for cat in cids.keys():
            cid_results[cat] = {}
            if self.sub_dirs is not None:
                cat_path = os.path.join(self.base_dir, self.sub_dirs[cat])

            for cid in cids[cat]:
                xml = self.template_query % {"id": cid,
                                             "query_type": "similarity",
                                             "type_specific": threshold_tag,
                                             "max_records": str(max_records)}
                page = urllib.request.urlopen(self.base_url, data=xml)
                output = page.read()
                page.close()
                filename = str(cid) + "_similarity."
                filepath = os.path.join(cat_path, filename)
                cid_results[cid] = self.check_and_process(output, filepath)

        return cid_results

    def scrape_super(self, cids, match_isotopes=False, match_charges=False,
                     match_tautomers=False,
                     rings_not_embedded=False,
                     single_double_bonds_match=True,
                     chains_match_rings=True,
                     strip_hydrogen=False,
                     stereo="ignore",
                     max_records=10000):
        """
        Generalized function for superstructure searches (searching for
        molecules that contain a given molecule within them) that either calls
        scrape_similar_rest or scrape_similar_pug, depending on if REST is being
        used. Kwargs are those used by PUG and PUG-REST for similarity queries.

        Parameter descriptions are largely taken from
        http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :param match_isotopes: Atoms must be of the same specified isotope.
        :param match_charges: Atoms must match the specified charge.
        :param match_tautomers: Allows matching with tautomers.
        :param rings_not_embedded: Rings may not be embedded in a larger system.
        :param single_double_bonds_match: In an aromatic compound, either single
        or double bonds may match the aromatic bonds.
        :param chains_match_rings: Chain bonds in the query may match rings in hits.
        :param strip_hydrogen: Remove explicit hydrogens before searching.
        :param stereo: How to handle stereoisomers: either "ignore", "exact",
        "relative", or "nonconflicting".
        :param max_records: Maximum number of hits.
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        if self.use_rest:
            self.scrape_super_rest(cids, match_isotopes=match_isotopes,
                                   match_charges=match_charges,
                                   match_tautomers=match_tautomers,
                                   rings_not_embedded=rings_not_embedded,
                                   single_double_bonds_match=single_double_bonds_match,
                                   chains_match_rings=chains_match_rings,
                                   strip_hydrogen=strip_hydrogen,
                                   stereo=stereo,
                                   max_records=max_records)
        else:
            self.scrape_super_pug(cids, match_isotopes=match_isotopes,
                                   match_charges=match_charges,
                                   match_tautomers=match_tautomers,
                                   rings_not_embedded=rings_not_embedded,
                                   single_double_bonds_match=single_double_bonds_match,
                                   chains_match_rings=chains_match_rings,
                                   strip_hydrogen=strip_hydrogen,
                                   stereo=stereo,
                                   max_records=max_records)

    def scrape_super_rest(self, cids, match_isotopes=False, match_charges=False,
                             match_tautomers=False,
                             rings_not_embedded=False,
                             single_double_bonds_match=True,
                             chains_match_rings=True,
                             strip_hydrogen=False,
                             stereo="ignore",
                             max_records=10000):
        """
        Searches by superstructure (searching for molecules that contain a given
        molecule within them) using PUG-REST.

        Parameter descriptions are largely taken from
        http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :param match_isotopes: Atoms must be of the same specified isotope.
        :param match_charges: Atoms must match the specified charge.
        :param match_tautomers: Allows matching with tautomers.
        :param rings_not_embedded: Rings may not be embedded in a larger system.
        :param single_double_bonds_match: In an aromatic compound, either single
        or double bonds may match the aromatic bonds.
        :param chains_match_rings: Chain bonds in the query may match rings in hits.
        :param strip_hydrogen: Remove explicit hydrogens before searching.
        :param stereo: How to handle stereoisomers: either "ignore", "exact",
        "relative", or "nonconflicting".
        :param max_records: Maximum number of hits.
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        output = {}
        queries_run = 0

        for cat in cids.keys():
            output[cat] = {}
            for cid in cids[cat]:
                queries_run = self.check_queries(queries_run)
                result = pcp.get_cids(cid, searchtype="superstructure",
                                      match_isotopes=match_isotopes,
                                      match_charges=match_charges,
                                      match_tautomers=match_tautomers,
                                      rings_not_embedded=rings_not_embedded,
                                      single_double_bonds_match=single_double_bonds_match,
                                      chains_match_rings=chains_match_rings,
                                      strip_hydrogen=strip_hydrogen,
                                      stereo=stereo,
                                      max_records=max_records, record_type="3d")
                output[cat][cid] = result

        return output

    def scrape_super_pug(self, cids, match_isotopes=False, match_charges=False,
                             match_tautomers=False,
                             rings_not_embedded=False,
                             single_double_bonds_match=True,
                             chains_match_rings=True,
                             strip_hydrogen=False,
                             stereo="ignore",
                             max_records=10000):
        pass