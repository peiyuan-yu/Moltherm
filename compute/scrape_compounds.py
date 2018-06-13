import os
import time
import re

# Unfortunately, not super useful, because of timeouts
import pubchempy as pcp

import urllib3
from bs4 import BeautifulSoup

import pybel

# Notes:
# record_type should be 3d, if possible (2d by default)
# Could search by superstructure (default options should be good, except maybe ChainsMatchRings - defaults true)
# Get SMILES string of cyclohexene, use that as superstructure
# Could search by similarity (play around with Threshold - default 90)

cids = {"products": [8079],
        "dienes": [7845],
        "dienophiles": [6325]
        }

diene_dir = "/Users/espottesmith/data/dienes/"
dienophile_dir = "/Users/espottesmith/data/dienophiles/"
product_dir = "/Users/espottesmith/data/products/"


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
                      %(type_specific_string)s
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
          <PCT-Download_format value="sdf"/>
          <PCT-Download_compression value="gzip"/>
        </PCT-Download>
      </PCT-InputData_download>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data>"""

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

    def download_files(self, cids):
        """

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :return:
        """

        pass

    def store_cids(self, cids):
        """

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :return:
        """

        pass

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
        pass

    def scrape_similar_pug(self, cids, threshold=90, max_records=10000):
        pass

    def scrape_super(self, cids, match_isotopes=False, match_charges=False,
                     match_tautomers=False,
                     rings_not_embedded=False,
                     single_double_bonds_match=True,
                     chains_match_rings=True,
                     strip_hydrogen=False,
                     stereo="ignore",
                     max_records=10000):
        """
        Generalized function that either calls scrape_similar_rest
        scrape_similar_pug, depending on if REST is being used. Kwargs are
        those used by PUG and PUG-REST for similarity queries.

        :param cids: A dict {"category": [ids]}, where each category is a
        molecular type of interest.
        :param match_isotopes:
        :param match_charges:
        :param match_tautomers:
        :param rings_not_embedded:
        :param single_double_bonds_match:
        :param chains_match_rings:
        :param strip_hydrogen:
        :param stereo:
        :param max_records:
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
        pass

    def scrape_super_pug(self, cids, match_isotopes=False, match_charges=False,
                             match_tautomers=False,
                             rings_not_embedded=False,
                             single_double_bonds_match=True,
                             chains_match_rings=True,
                             strip_hydrogen=False,
                             stereo="ignore",
                             max_records=10000):
        pass