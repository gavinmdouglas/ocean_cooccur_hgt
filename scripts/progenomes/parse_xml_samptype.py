#!/usr/bin/python3
import os
import glob
import xml.etree.ElementTree as ET
import sys

attributes_of_interest = ['sample-type', 'isolation-source', 'isolate', 'completeness_estimated', 'contamination_estimated', 'completeness', 'contamination', 'completeness-score', 'contamination-score', 'metagenome-source', 'metagenomic']

def parse_biosample_types(xml_file):
    """
    Parse NCBI Biosample XML file to extract sample types.
    
    Args:
        xml_file (str): Path to XML file containing biosample data
        
    Returns:
        dict: Mapping of biosample accessions to their sample types
    """
    
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Dictionary to store results
    sample_info = {}
    
    # Iterate through each BioSample
    for biosample in root.findall(".//BioSample"):
        # Get the accession
        accession = biosample.get('accession')

        sample_info[accession] = {}
        
        # Look for sample-type attribute
        
        for attribute in attributes_of_interest:
            sample_info[accession][attribute] = "NA"

        attributes = biosample.find('.//Attributes')
        if attributes is not None:
            for attribute in attributes.findall('Attribute'):
                attribute_name = attribute.get('attribute_name')
                attribute_name = attribute_name.replace(' ', '-')
                if attribute_name in attributes_of_interest:
                    sample_info[accession][attribute_name] = attribute.text.replace(' ', '_').replace('\t', '_')
    
        # Get "Status" as well.
        status_elem = biosample.find('.//Status')
        status = status_elem.get('status')
        if status:
            sample_info[accession]['status'] = status
        else:
            sample_info[accession]['status'] = 'NA'
    return sample_info


biosample_ids = set()
with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/all_biosample_ids.txt', 'r') as id_fh:
    for line in id_fh:
        biosample_ids.add(line.strip())

with open('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/additional_biosample_map.txt', 'r') as id_fh:
    for line in id_fh:
        line = line.strip().split()
        biosample_ids.add(line[1])

xml_files = glob.glob(os.path.join('/mfs/gdouglas/projects/ocean_mags/progenomes_analyses/genome_prep/xml_sra_metadata', "*.xml"))

all_parsed = {}
for xml_file in xml_files:
    all_parsed.update(parse_biosample_types(xml_file))

header = ['biosample_id', 'Presence', 'Status']
header.extend(attributes_of_interest)
print('\t'.join(header))

for biosample_id in biosample_ids:
    outline = [biosample_id]
    if biosample_id in all_parsed:
        outline.append('Present')
        outline.append(all_parsed[biosample_id]['status'])
        for attribute in attributes_of_interest:
            outline.append(all_parsed[biosample_id][attribute])
    else:
        outline.append('Missing')
        outline.extend(['NA'] * (len(attributes_of_interest) + 1))
    
    print('\t'.join(outline))
