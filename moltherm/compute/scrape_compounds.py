# coding: utf-8

import os
import time
from io import StringIO
import urllib.request
import requests

import pubchempy as pcp
from chemspipy import ChemSpider
from chemspipy.errors import *

from bs4 import BeautifulSoup

from monty.json import jsanitize
from pymongo import MongoClient

from matgendb.dbconfig import DBConfig
from pymatgen.io.babel import BabelMolAdaptor

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "June 2018"


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
                  <PCT-QueryCompoundCS_results>300</PCT-QueryCompoundCS_results>
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
                os.makedirs(self.base_dir)
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

    def check_and_process(self, output, request_type, path=None):
        """
        Parses output XML. If the output returns a file, downloads that file
        at the specified path; if the output returns a request id, continue to
        check at that id until the request is complete.

        :param output: An XML file generated by download_files_pug
        :param request_type: Either "download" or "query"
        :param path: A path string for downloads, generated by
        download_files_pug
        :return:
        """

        wait_time = 1
        num_checks = 500

        first_check = True

        for i in range(num_checks):
            if request_type == "download":
                parsed = BeautifulSoup(output, "lxml-xml")
                finished = parsed.find("PCT-Download-URL_url")

                if finished is not None:
                    ftp_url = finished.text
                    opener = urllib.request.URLopener()
                    opener.retrieve(ftp_url, path)
                    return

                request_id = parsed.find("PCT-Request_reqid").text

            elif request_type == "query":

                if first_check:
                    parsed = BeautifulSoup(output, "html.parser")
                    request_id = parsed.find("ListKey").text
                    first_check = False
                else:

                    parsed = BeautifulSoup(output, "lxml-xml")
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

                    request_id = parsed.find("PCT-Request_reqid").text

            time.sleep(wait_time)

            request = self.template_check % {"request_id": request_id}
            request_file = StringIO()
            request_file.write(request)

            page = urllib.request.urlopen(self.base_url, data=request_file)
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

        order = 0

        formats = ["SDF"]
        if pngs:
            formats.append("PNG")

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

            for format in formats:
                for cid in download_ids:
                    filename = str(cid) + "_" + str(order) + "." + format.lower()
                    filepath = os.path.join(cat_path, filename)
                    pcp.download(format, filepath, cid, overwrite=True)
                    order += 1

    def download_files_pug(self, cids, pngs=True, download_parents=False):
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
                xml_file = StringIO()
                xml_file.write(xml)
                output = urllib.request.urlopen(self.base_url, data=xml_file)
                filename = str(cat) + "_." + format + ".gz"
                filepath = os.path.join(cat_path, filename)
                self.check_and_process(output, "download", path=filepath)

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
            return self.scrape_similar_rest(cids, threshold=threshold,
                                     max_records=max_records)
        else:
            return self.scrape_similar_pug(cids, threshold=threshold,
                                    max_records=max_records)

    def scrape_similar_rest(self, cids, threshold=90, max_records=10000):
        """
        Searches by similarity (2D) using PUG-REST.

        :param cids:
        :param threshold:
        :param max_records:
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        output = {}
        queries_run = 0

        for cat in cids.keys():
            output[cat] = {}
            for cid in cids[cat]:
                queries_run = self.check_queries(queries_run)
                result = pcp.get_cids(cid, namespace="cid",
                                           domain="compound",
                                           searchtype="similarity",
                                           threshold=threshold,
                                           max_records=max_records)
                output[cat][cid] = result

        return output

    def scrape_similar_pug(self, cids, threshold=90, max_records=10000):
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

        cid_results = {}

        threshold_tag = """"<PCT-CSSimilarity_threshold>
        %(threshold)s
        </PCT-CSSimilarity_threshold>""" % {"threshold": str(threshold)}

        for cat in cids.keys():
            cid_results[cat] = {}

            for cid in cids[cat]:
                xml = self.template_query % {"id": cid,
                                             "query_type": "similarity",
                                             "type_specific": threshold_tag,
                                             "max_records": str(max_records)}
                xml_file = StringIO()
                xml_file.write(xml)
                page = urllib.request.urlopen(self.base_url, data=xml_file)
                output = page.read()
                page.close()

                cid_results[cid] = self.check_and_process(output, "query")

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
            return self.scrape_super_rest(cids, match_isotopes=match_isotopes,
                                   match_charges=match_charges,
                                   match_tautomers=match_tautomers,
                                   rings_not_embedded=rings_not_embedded,
                                   single_double_bonds_match=single_double_bonds_match,
                                   chains_match_rings=chains_match_rings,
                                   strip_hydrogen=strip_hydrogen,
                                   stereo=stereo)
        else:
            return self.scrape_super_pug(cids)

    def scrape_super_rest(self, cids, match_isotopes=False, match_charges=False,
                          match_tautomers=False,
                          rings_not_embedded=False,
                          single_double_bonds_match=True,
                          chains_match_rings=True,
                          strip_hydrogen=False,
                          stereo="ignore",
                          max_records=10000):
        """
        Generalized function for superstructure searches (searching for
        molecules that contain a given molecule within them).
        Kwargs are those used by PUG and PUG-REST for similarity queries.

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
                result = pcp.get_cids(cid, namespace="cid",
                                      domain="compound",
                                      searchtype="superstructure",
                                      match_isotopes=match_isotopes,
                                      match_charges=match_charges,
                                      match_tautomers=match_tautomers,
                                      rings_not_embedded=rings_not_embedded,
                                      single_double_bonds_match=single_double_bonds_match,
                                      chains_match_rings=chains_match_rings,
                                      strip_hydrogen=strip_hydrogen,
                                      stereo=stereo,
                                      max_records=max_records)
                output[cat][cid] = result

        return output

    def scrape_super_pug(self, cids):
        """
        Generalized function for superstructure searches (searching for
        molecules that contain a given molecule within them).

        NOTE: I cannot find documentation on this kind of search online. As such,
        I'm only even attempting to create the most basic possible XML file for
        this type of search.
        TODO: Figure out what the schema is actually supposed to look like.

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :return: A dict {"category": {id:[matches]}}, where each id in each
        category is stored along with its matches (which take the form of CIDs.
        """

        cid_results = {}

        for cat in cids.keys():
            cid_results[cat] = {}

            for cid in cids[cat]:
                xml = self.template_query % {"id": cid,
                                             "query_type": "superstructure",
                                             "type_specific": "",
                                             "max_records": 1000000}
                xml_file = StringIO()
                xml_file.write(xml)
                page = urllib.request.urlopen(self.base_url, data=xml_file)
                output = page.read()
                page.close()

                cid_results[cid] = self.check_and_process(output, "query")

        return cid_results


class ReaxysScraper:
    """
    Not so much a scraper as a parser (for now). Takes results from Reaxys
    query. Can convert to Reaction namedtuple, which includes important
    metadata and CTAB data. Can also turn that CTAB data into CTAB files. This
    allows easy interoperability with OpenBabel and pymatgen.
    """

    def __init__(self, base_dir):
        """
        ReaxysScraper

        :param base_dir: Base directory for parsing and generating data.
        """

        self.base_dir = base_dir

    def parse_reaxys_xml(self, filename):
        """
        Parses an XML file generated by the Reaxys API.

        :param filename: str referring to XML file from Reaxys.
        :return: List of dicts including reactant CTAB, product CTAB, and
        metadata.
        """

        results = []

        filepath = os.path.join(self.base_dir, filename)

        with open(filepath, 'r') as fileobj:
            xml = fileobj.read()
            parsed = BeautifulSoup(xml, "lxml-xml")

            reactions = parsed.find_all("reaction")

            for reaction in reactions:
                # Screen for reactions with more than two reactants
                # or more than one product
                if len(reaction.find_all("RY.PRO")) > 2\
                        or len(reaction.find_all("RY.RCT")) > 2:
                    continue

                # Generate metadata from reaction header information
                # Will be passed along with CTAB information
                index = int(reaction["index"])
                rxn_id = int(reaction.find("RX.ID").text)
                solvents = set([sol.text for sol
                                in reaction.find_all("RXD.SOL")])

                rct_ids = reaction.find_all("RX.RXRN")
                rct_names = reaction.find_all("RX.RCT")
                rct_meta = [(int(e.text), rct_names[i].text) for i, e
                            in enumerate(rct_ids)]

                pro_ids = reaction.find_all("RX.PXRN")
                pro_names = reaction.find_all("RX.PRO")
                pro_meta = [(int(e.text), pro_names[i].text) for i, e
                            in enumerate(pro_ids)]

                meta = {"index": index,
                        "rxn_id": rxn_id,
                        "solvents": solvents,
                        "rct_meta": sorted(rct_meta, key=lambda x: x[0]),
                        "pro_meta": sorted(pro_meta, key=lambda x: x[0])}

                # Capture reactant CTAB information
                # Make sure that ordering is the same for metadata and CTAB
                rcts = sorted(reaction.find_all("RY.RCT"),
                              key=lambda x: int(x["rn"]))
                rcts = [rct.text for rct in rcts]

                pros = sorted(reaction.find_all("RY.PRO"),
                              key=lambda x: int(x["rn"]))

                pros = [pro.text for pro in pros]

                rxn = {"rcts": rcts, "pros": pros, "meta": meta}

                results.append(rxn)

        return results

    def store_reaxys_reactions_files(self, reactions):
        """
        Create CTAB (.mol) files from parsed XML data.

        :param reactions: List of reactions, defined as above.
        :return:
        """

        for reaction in reactions:

            # Label directories by index and reaction id
            dir_name = str(reaction["meta"]["index"]) + "_" + \
                       str(reaction["meta"]["rxn_id"])
            path = os.path.join(self.base_dir, "reactions", dir_name)

            if not os.path.exists(path):
                os.makedirs(path)

            # Create metadata file
            with open(os.path.join(path, "meta.xml"), 'w') as file:
                file.write("<metadata>")
                file.write("<index>%(index)s</index>\n" % {"index": str(reaction["meta"]["index"])})
                file.write("<reaxysid>%(id)s</reaxysid>\n" % {"id": str(reaction["meta"]["rxn_id"])})

                file.write("<solvents>%(solvents)s</solvents>\n" % {"solvents": ",".join(reaction["meta"]["solvents"])})

                for i, rct in enumerate(reaction["meta"]["rct_meta"]):
                    rct_info = {"num": str(i),
                                "name": rct[1],
                                "id": str(rct[0])
                                }
                    file.write("""<rct num=%(num)s>
<rctname>%(name)s</rctname>
<rctid>%(id)s</rctid>
</rct>\n""" % rct_info)

                for i, pro in enumerate(reaction["meta"]["pro_meta"]):
                    pro_info = {"num": str(i),
                                "name": pro[1],
                                "id": str(pro[0])}
                    file.write("""<pro>
<proname>%(name)s</proname>
<proid>%(id)s</proid>
</pro>\n""" % pro_info)

                file.write("</metadata>")

            # Create reactant files, named with their Reaxys IDs
            reactants = reaction["meta"]["rct_meta"]
            for i, e in enumerate(reactants):
                filename = "rct_" + str(i) + "_" + str(e[0]) + ".mol"
                with open(os.path.join(path, filename), 'w') as file:
                    file.write(reaction["rcts"][i])

            # Create product file, named with its Reaxys ID
            products = reaction["meta"]["pro_meta"]
            for i, e in enumerate(products):
                filename = "pro_" + str(i) + ".mol"
                with open(os.path.join(path, filename), 'w') as file:
                    file.write(reaction["pros"][i])

    @staticmethod
    def store_reaxys_reactions_db(reactions, dbfile="db.json",
                                  collection="reactions"):
        """
        Insert reaction information into a MongoDB database.

        :param reactions: List of reactions, defined as above.
        :param dbfile: A config file indicating the database into which the
        reactions will be inserted.
        :param collection: Collection within the database in which to store
        the reactions.
        :return: List of rxn_ids that were
        """

        # Set up MongoDB database with pymatgen-db and pymongo
        config = DBConfig(config_file=dbfile)

        client = MongoClient(config.host, config.port, username=config.user,
                             password=config.password)

        db = client[config.dbname]
        collection = db[collection]

        just_added = []

        if "rxn_id" not in collection.index_information():
            collection.create_index("rxn_id", unique=True)

        for reaction in reactions:
            try:
                reaction["rxn_id"] = reaction["meta"]["rxn_id"]
                collection.insert_one(reaction)
                just_added.append(reaction["rxn_id"])
            except:
                continue

        return just_added


class ChemSpiderScraper:
    """
    Uses ChemSpiPy as well as BeautifulSoup to extract boiling and melting
    points for chemicals of interest.
    """

    def __init__(self, token, base_dir, subdirs=False):
        """
        Initialize ChemSpider Scraper.

        :param token: ChemSpider API token (str). This should be obtained from
        ChemSpider ()
        :param base_dir: Directory where information should be stored.
        :param subdirs: Should subdirectories be used? By default, this is False,
        meaning that all files will be stored in the base_dir
        """

        self.token = token
        self.base_dir = base_dir
        self.subdirs = subdirs

        self.chemspi = ChemSpider(self.token)

        # For cases where we cannot use ChemSpiPy, we need to be able to set
        # up requests
        self.base_url = "https://api.rsc.org/compounds/v1/"
        self.headers = {"apikey": self.token, "Content-Type": "application/json"}

    def get_chemspider_ids(self, molecules, max_attempts=100):
        """
        From list of :class: pymatgen.core.structure.Molecule, extract
        ChemSpider IDs for use in subsequent searches.

        :param molecules: List of Molecules.
        :param max_attempts: How many times should a query be checked before
            it should be given up?
        :return: dict {molecule: ids}, where ids is an list of ints.
        """

        results = {}

        for molecule in molecules:
            adaptor = BabelMolAdaptor(molecule)
            # Use Canonical SMILES to ensure uniqueness
            smiles = adaptor.pybel_mol.write("can").strip()

            try:
                # Attempt to use the ChemSpiPy package
                # This significantly simplifies the work to be done
                results[smiles] = self.chemspi.simple_search(smiles)
            except ChemSpiPyServerError:
                init_url = self.base_url + "filter/smiles"
                data = {"smiles": str(smiles)}
                init_req = requests.post(init_url, json=data, headers=self.headers)
                try:
                    query_id = init_req.json()['queryId']
                except KeyError:
                    raise RuntimeError("Response did not include queryId key!")

                status_url = self.base_url + "filter/{}/status".format(query_id)
                for i in range(max_attempts):
                    status_request = requests.get(status_url, headers=self.headers)

                    if status_request.json()['status'].lower() == 'complete':
                        results_url = self.base_url + \
                                      "filter/{}/results".format(query_id)

                        results_request = requests.get(results_url,
                                                       headers=self.headers)

                        results[smiles] = results_request.json().get("results", [])
                        break

        return results

    def extract_boiling_point(self, csids):
        """
        Scrape experimental and/or predicted boiling point information from
        ChemSpider.

        :param csids: list of ChemSpider id.
        :return: dict containing all listed boiling points
        """

        results = {}

        for csid in csids:
            results[csid] = {}

            url = """http://parts.chemspider.com/JSON.ashx?op=GetRecordsAsCompounds&csids[0]={}&serfilter=Compound[PredictedProperties|ExperimentalProperties]""".format(str(csid))
            request = requests.get(url)

            parsed = BeautifulSoup(request.text, "lxml")

            # Turn text into dict of lists of dicts
            data = eval(parsed.body.p.text)[0]

            results[csid]["pred"] = []
            for pred in data["PredictedProperties"]:
                if pred["Name"] == "Boiling Point":
                    results[csid]["pred"].append({
                        "units": pred.get("Units", None),
                        "value": pred.get("Value", None)})

            results[csid]["exp"] = []
            for exp in data["ExperimentalProperties"]:
                if exp["Name"] == "Experimental Boiling Point":
                    results[csid]["exp"].append({
                        "units": exp.get("Units", None),
                        "value": exp.get("Value", None),
                        "source": exp.get("DataSourceName", None)})

        return results


    def extract_melting_point(self, csids):
        """
        Scrape experimental and/or predicted boiling point information from
        ChemSpider.

        :param csids: list of ChemSpider ids.
        :return: dict containing all listed boiling points
        """

        results = {}

        for csid in csids:
            results[csid] = {}

            url = """http://parts.chemspider.com/JSON.ashx?op=GetRecordsAsCompounds&csids[0]={}&serfilter=Compound[PredictedProperties|ExperimentalProperties]""".format(str(csid))
            request = requests.get(url)

            parsed = BeautifulSoup(request.text, "lxml")

            # Turn text into dict of lists of dicts
            data = eval(parsed.body.p.text)[0]

            results[csid]["pred"] = []
            for pred in data["PredictedProperties"]:
                if pred["Name"] == "Boiling Point":
                    results[csid]["pred"].append({
                        "units": pred.get("Units"),
                        "value": pred.get("Value")})

            results[csid]["exp"] = []
            for exp in data["ExperimentalProperties"]:
                if exp["Name"] == "Experimental Melting Point":
                    results[csid]["exp"].append({
                        "units": exp.get("Units"),
                        "value": exp.get("Value"),
                        "source": exp.get("DataSourceName")})

        return results

