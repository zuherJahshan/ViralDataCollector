import logging
import os

from typing import List, Dict
Accession = str
Lineage = str
Sequence = str
PropertyName = str
Value = str
AccessionsPropertiesMap = Dict[Accession, Dict[PropertyName, Value]]

class Covid19DataPortalAccessionFetcher(object):
    def __init__(self,
                 data_path = "data/raw/",
                 log_file_name = "file.log"):

        # Define the local parameters
        self.log_file_name = log_file_name
        self.data_path = data_path
        if self.data_path[-1] != '/':
            self.data_path += '/'

        os.makedirs(self.data_path,
                    exist_ok=True)

        # Check for the accessions file, if not found return an error
        if not os.path.exists("accessions.tsv"):
            logging.critical("the \"accessions.tsv\" file does not exist.\n")
            exit(-1)

        # build the AccessionsPropertiesMap
        self.accessions_properties_map = self._buildAccessionsPropertiesMap()

        # Build the LineageAccessionsMap
        self.lineage_accessions_map = self._buildLineageAccessionsMap()

        # Build the local accessions set
        if os.path.exists(self.data_path):
            self.local_accessions_set = set([file.split('.')[0] for file in os.listdir("data/raw")])
        else:
            self.local_accessions_set = set()

    # Interface
    def downloadAccessions(self,
                           accessions: List[Accession],
                           download_all=True):
        all_accessions_exist = True
        existing_accessions = []
        for accession in accessions:
            if accession not in self.accessions_properties_map:
                all_accessions_exist = False
                print("accession {} does not exist in the system.".format(accession))
            elif accession in self.local_accessions_set:
                continue
            else:
                existing_accessions.append(accession)
        if not all_accessions_exist and download_all:
            print("Accessions will not be downloaded. To download only existing accessions,\
            please specify download_all=False when calling this function.")
            return

        valid_download = True
        for accession in accessions:
            self._downloadAccession(accession)
            valid_download &= self._checkDownload(accession)
        if not valid_download:
            print("not all accessions has successfully been downloaded, for more information check the log.")

    def getLineagesList(self):
        lineages = []
        for lineage in self.lineage_accessions_map:
            lineages.append(lineage)
        return lineages

    def getAccessionsByLineage(self,
                               lineage: Lineage) -> List[Accession]:
        if not lineage in self.lineage_accessions_map:
            print("Lineage {} is not specified".format(lineage))
            return []
        return list(self.lineage_accessions_map[lineage])

    ################################
    ####### Private methods ########
    ################################
    def _buildAccessionsPropertiesMap(self):
        accessions_file = open("accessions.tsv")
        properties_line = True
        accessions_properties_map = {}
        properties = []
        for line in accessions_file:
            def getElementsInLine(line):
                elements = line[:-1].split('\t')
                return [element.replace("\"", "") for element in elements]
            if properties_line:
                # extract the properties names
                properties = getElementsInLine(line)[1:] # The [1:] is for the sake of discarding the acc property
                properties_line = False
            else:
                # Extract from the line the properties values and the accession
                elements_in_line = getElementsInLine(line)
                accession = elements_in_line[0]
                properties_values = elements_in_line[1:]

                # Build the Info map according to the different properties
                properties_values_map = {}
                for i in range(len(properties)):
                    properties_values_map.update({properties[i]: properties_values[i]})

                # Append to the AccessionsPropertiesMap the pair accessions properties
                accessions_properties_map.update({accession: properties_values_map})
        accessions_file.close()
        return accessions_properties_map

    def _buildLineageAccessionsMap(self):
        lineage_accession_map = {}
        # Loop over the accessions
        for accession in self.accessions_properties_map:
            # Find the lineage of the accession
            lineage = self.accessions_properties_map[accession]["lineage"]

            # If lineage exists in the map, insert the accession to the existing list, else create a new entry in the map.
            if lineage in lineage_accession_map:
                lineage_accession_map[lineage].update({accession})
            else:
                lineage_accession_map.update({lineage: {accession}})
        return lineage_accession_map


    def _downloadAccession(self,
                           accession: Accession):
        os.system("wget https://www.ebi.ac.uk/ena/browser/api/fasta/{}.1?download=true -O {}{}.fasta".format(accession,
                                                                                                             self.data_path,
                                                                                                             accession))

    def _checkDownload(self,
                      accession: Accession) -> bool:
        if not os.path.exists("data/{}.fasta".format(accession)):
            logging.warning("accession {} found in the accessions.tsv but could not be downloaded.\n")
            return False
        self.local_accessions_set.add(accession)
        return True
