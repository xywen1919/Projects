#!/usr/bin/env python
# coding: utf-8

# In[3]:


# import packages
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Seq import Seq
from bioservices import *
from IPython.display import Image, HTML  # Show images inline
import shutil  # Standard library packages
import io  # Standard library packages
import os  # Standard library packages
import glob  # Standard library packages
import shutil  # Standard library packages
from os import path
import re  # for regular expression operations
import pandas as pd  # for usage of dataframe

from Bio.Graphics.ColorSpiral import ColorSpiral
import random


# # Define functions

# In[4]:


# define working environment
file_path = "/Users/wen_x/Downloads/BioData/FinalProject"
if not os.path.exists(file_path):  # create a project directory folder
    os.makedirs(file_path)

# input email in order to use entrez
my_email = "xw2470@nyu"
Entrez.email = my_email

# created KEGG object
k = KEGG(verbose=False)


# In[3]:


# define function to fetch disease entry number given disease name
def get_diseaseEntryNumber(diseaseName):
    """
    input: diseaseName 
	return: disease entry_number
    """
    diseaseEntry = k.find(database='disease', query=diseaseName).split('\t')[0]
    entry_number = diseaseEntry.split(':')[1]
    return entry_number


# In[4]:


# function to write annotation from entry number to local folder 
def get_annotation(entry_number):
    """
    input: entry_number 
	return: the corresponding annotation and write to local folder
    """
    # read content for each entry_number & write to
    target_fileName = file_path + "/annotation_%s.txt" % entry_number  # file name
    if not os.path.isfile(target_fileName):  # create file if not exists
        open(target_fileName, 'w').write(k.get(entry_number))  # write to file
        print("The annotation of %s is saved." % entry_number)  # print message


# In[5]:


# function to fetch pathway ID from disease entry number
def get_pathwayID(disease_entry_number):
    """
    input: disease entry number
	return: list of pathway ID     
    """
    # get annotation from entry_number
    get_annotation(disease_entry_number)
    # read input file & find pathway_ID
    input_file = file_path + "/annotation_%s.txt" % disease_entry_number
    with open(input_file, 'r') as f:
        pattern = re.compile(r"(hsa\d{5})")
        pathwayID = re.findall(pattern, f.read())
    # return value
    return pathwayID


# In[6]:


# function to write the pathway_map file to local folder 
def draw_kegg_map(disease_entry_number):
    """
    input: disease entry number
	return: render a local PDF file of a KEGG map 
    """
    # get pathwayID from entry_number
    pathway_id_list = get_pathwayID(disease_entry_number)
    for pathway_id in pathway_id_list:
        # create fileName for image saving
        img_filename = file_path + "/pathwayMap_%s_%s.pdf" % (disease_entry_number, pathway_id)
        if not os.path.isfile(img_filename):  # check existence of file
            pathway = KGML_parser.read(k.get(pathway_id, "kgml"))  # Get the map info
            canvas = KGMLCanvas(pathway)  # draw pathway map
            canvas.import_imagemap = True  # add background map
            canvas.draw(img_filename)  # write to target img_file
            # print out message
            print("Pathway map of %s_%s has been saved." % (disease_entry_number, pathway_id))


# In[7]:


# function to fetch the network ID from disease entry ID
def get_networkID(disease_entry_number):
    """
    input: disease entry_number 			
    return: list of network id
    """
    # get annotation from disease_entry_number
    get_annotation(disease_entry_number)
    # read input file & find network_IDs
    input_file = file_path + "/annotation_%s.txt" % disease_entry_number
    with open(input_file, 'r') as f:
        pattern = re.compile(r"(nt\d{5})")
        networkID = re.findall(pattern, f.read())
    # return values
    return networkID


# In[8]:


# function to write the annotation of networks to local folder 
def get_networkAnnotation(disease_entry_number):
    """
    input: disease entry_number 			
    return: network annotation .txt file
    """
    for i in get_networkID(disease_entry_number):
        network_entry = "network:%s" % i
        target_fileName = file_path + "/annotation_%s.txt" % i  # file name
        if not os.path.isfile(target_fileName):  # create file if not exists
            open(target_fileName, 'w').write(k.get(network_entry))  # write to file
            print("The annotation of %s is saved." % network_entry)  # print message


# In[9]:


# function to fetch element ID from disease entry number
def get_elementID(disease_entry_number):
    """
    input: disease entry_number -> annotations
			-> network id annotations			
    return: list of element id 
    """
    # write disease annotation to local folder
    get_annotation(disease_entry_number)
    # read disease input file from local folder & find element_IDs
    input_file = file_path + "/annotation_%s.txt" % disease_entry_number
    with open(input_file) as f:
        pattern = re.compile(r"(N\d{5})")
        elementID_diseAnnot = re.findall(pattern, f.read())

    # write network annotation to local folder
    get_networkAnnotation(disease_entry_number)
    # read network input file from local folder & find element_IDs
    elementID_ntwkAnnot = []
    for input_file in glob.glob(file_path + "annotation_nt*.txt"):
        with open(input_file) as f:
            pattern = re.compile(r"(N\d{5})")
            elementID_ntwkAnnot.extend(re.findall(pattern, f.read()))

    # combine elementID list and remove duplicates
    elementID = list(dict.fromkeys(elementID_diseAnnot + elementID_ntwkAnnot))
    # return value
    return elementID


# In[10]:


# function to write list of geneIDs to local folder 
def get_geneID(disease_entry_number):
    """
    input: disease entry_number 
			-> list of element_ID & pathway_ID 
	return: list of geneID & .txt file
    """
    # get list of pathwayID + elementID from disease_entry_number
    input_id_list = get_pathwayID(disease_entry_number) + get_elementID(disease_entry_number)
    # loop read content for each id, extract gene id and add to list
    pattern = re.compile(r"\s(\d{3,9})\s")
    gene_id = []
    for i in input_id_list:
        gene_id.extend(re.findall(pattern, k.get(i)))
    geneID = list(dict.fromkeys(gene_id))  # Remove duplicates
    # define target file name; create new file if not exists
    target_fileName = file_path + "/geneID_%s.txt" % disease_entry_number
    if not os.path.isfile(target_fileName):  # check if the file exists
        # write to file
        open(target_fileName, 'w').write(str(geneID))
    # return value
    return geneID


# In[11]:


# function to fetch metabolite id 
def get_metaboliteID(disease_entry_number):
    """
    input: disease_entry_number 
			-> elementID 
    return: list of metaboliteIDs
    """
    # get list of elementID from disease entry_number
    element_id = get_elementID(disease_entry_number)
    # define pattern and find all matches in each element id annotation
    pattern = re.compile(r"(C\d{5})")
    metabolite_id = []
    for i in element_id:
        metabolite_id.extend(re.findall(pattern, k.get(i)))
    metaboliteID = list(dict.fromkeys(metabolite_id))  # remove duplicates
    return metaboliteID


# In[12]:


# function to fetch perturbant id 
def get_perturbantID(disease_entry_number):
    """
    input: disease_entry_number 
			-> elementID  
    return: list of perturbantID
    """
    # get list of elementID from entry_number
    element_id = get_elementID(disease_entry_number)
    # define pattern and find all matches in each element id annotation
    pattern = re.compile(r"(K\d{5})")
    perturbant_id = []
    for i in element_id:
        perturbant_id.extend(re.findall(pattern, k.get(i)))
    perturbantID = list(dict.fromkeys(perturbant_id))  # remove duplicates
    return perturbantID


# In[14]:


# Function to fetch pubmed id 
def get_pubmedID(disease_entry_number):
   """
   input: disease entry number
		-> annotations of [disease_entry_number, pathwayID, elementID, perturbantID, metaboliteID]
   return: list of pubmedIDs 
   """
   # get input_id_list from disease_entry_number
   input_id_list = [str(disease_entry_number)]
   input_id_list.extend(get_pathwayID(disease_entry_number) + get_elementID(disease_entry_number))
   input_id_list.extend(get_perturbantID(disease_entry_number) + get_metaboliteID(disease_entry_number))
   # define pattern and find all matches in each id annotation
   pattern = re.compile(r"PMID:(\d*)\s")
   pubmed_id = []
   for i in input_id_list:
       pubmed_id.extend(re.findall(pattern, k.get(i)))
   pubmedID = list(dict.fromkeys(pubmed_id))  # remove duplicates
   return pubmedID


# In[35]:


# function to write disease reference to local folder
def get_reference(disease_entry_number):
    """
    input: disease_entry_number 
			-> pubmedID
    return: render a .txt file with pubmedID & title 
    """
    # define target file name; create new file if not exists
    target_fileName = file_path + "/reference_title_%s.txt" % disease_entry_number
    if not os.path.isfile(target_fileName):
        # write to file
        with open(target_fileName, 'a') as f:
            f.write('pubmedID\tTitle\n')  # column names
            # get list of pubmed id from disease_entry_number
            pubmed_id = get_pubmedID(disease_entry_number)
            for i in pubmed_id:  # loop to read pubmed id and title for each reference and write to file
                handle = Entrez.esummary(db='pubmed', id=i)
                record = Entrez.read(handle)
                f.write(record[0]['Id'] + '\t' + record[0]['Title'] + '\n')
                handle.close()
            


# In[16]:


# function to get 'fasta' files from designated database and keyword
def get_fasta_database(databaseName, keyword):
    """
    input: identified databaseName & keyword
    return: render a .txt file with 200 fasta seq
    """
    # create fasta file
    db_record = Entrez.esearch(db=databaseName, term=keyword, retmax=200)
    id_list = Entrez.read(db_record)["IdList"]
    db_fasta = Entrez.efetch(db=databaseName, id=id_list, rettype='fasta', retmode='txt')
    # create file if not exists
    target_fileName = file_path + "/%s_%s_fasta.txt" % (keyword, databaseName)
    if not os.path.isfile(target_fileName):
        # write to file
        open(target_fileName, 'w').write(db_fasta.read())
        # printout message after complete writing
        print("The %s %s database has been save." % (keyword, databaseName))
    return db_fasta


# In[5]:


def summarize_fasta(*files):
    """    
    :param files: 
    :return: id dictionary & print out record.description
    """	
    id_dict = {}   # define a dictionary to hold id list
    for x in files:  # loop fasta files to record id and print out record.description
        id_dict[x] = []
        records = SeqIO.parse(x, 'fasta')  # define SeqIO objects
        for record in records:   # loop to read each record and print out items		
            id_dict[x].append(record.id)
            print(record.description)
        print("=====================")
    return id_dict


# In[17]:


# define a function to move files
def move_files(from_path, to_path, file_name):
    """
    :param from_path: 
    :param to_path: 
    :param file_name: 
    :return: file list in target folder
    """
    if not os.path.isfile(to_path + file_name):
        for x in glob.glob(from_path + file_name):
            shutil.move(x, to_path)
    filelist_in_folder = os.listdir(to_path)
    print("The files have been moved to %s." % to_path)
    return filelist_in_folder


# In[4]:


# define function to parse .txt blast results
def analyze_blast_result2(txt_file):
    """
    :param txt_file:
    :return: blast_parsed.txt in the local folder
    """
    infile = open(txt_file)  # open file to analyze
    target_file = file_path + '/blast_parsed.txt'  # target file name
    pattern = re.compile(r"(=|MZ)")
    with open(target_file, 'a') as f:
        for line in infile:
            if re.search(pattern, line):  # search for pattern line by line
                f.write(line)  # write to file if not match
    infile.close()  # close input analyze file
    print("Analysis of the Blast result has been written to local folder.")  # print out message


# # Analysis

# In[19]:


# get disease entry number for "COVID-19"
print(get_diseaseEntryNumber('COVID-19'))


# In[20]:


# Q1. write KEGG annotation file for "COVID-19" 
get_annotation(get_diseaseEntryNumber('COVID-19'))


# In[21]:


# gether list of pathwayID for "COVID-19"
print(get_pathwayID(get_diseaseEntryNumber('COVID-19')))


# In[22]:


# Q2. write pathway image .pdf file for "COVID-19"
draw_kegg_map(get_diseaseEntryNumber('COVID-19'))


# In[23]:


# gether list of networkID
print(get_networkID(get_diseaseEntryNumber('COVID-19')))


# In[24]:


# Q3. write annotation of networks as .txt file for "COVID-19"
print(get_networkAnnotation(get_diseaseEntryNumber('COVID-19')))


# In[25]:


# gether list of elementID
print(get_elementID(get_diseaseEntryNumber('COVID-19')))


# In[26]:


# Q4. return list of gene id and write to .txt file
print(get_geneID(get_diseaseEntryNumber('COVID-19')))


# In[27]:


# gether list of metabolite id
print(get_metaboliteID(get_diseaseEntryNumber('COVID-19')))


# In[28]:


# gether list of perturbant id
print(get_perturbantID(get_diseaseEntryNumber('COVID-19')))


# In[29]:


# gether pubmed id from disease entry number
print(get_pubmedID(get_diseaseEntryNumber('COVID-19')))


# In[38]:


#Q7. get pubmedID with its reference title written to .txt file
print(get_reference(get_diseaseEntryNumber('COVID-19')))


# In[39]:


# Q5. to get 200 genomic sequence for SARS-CoV-2 as blast database
print(get_fasta_database(databaseName='nucleotide', keyword='SARS-CoV-2'))


# In[40]:


# to get 200 NTsequence of S protein
print(get_fasta_database(databaseName='nucleotide', keyword='SARS-CoV-2 S protein'))


# In[6]:


# analyze sequence description
summarize_fasta(file_path+'/SARS-CoV-2 S protein_nucleotide_fasta.txt', file_path+'/SARS-CoV-2_nucleotide_fasta.txt')


# In[41]:


# move fasta files to cygwin64 home folder for blast
move_files(from_path=file_path, to_path='C:/cygwin64/home/wen_x', file_name='/SARS-CoV-2 S protein_nucleotide_fasta.txt')
move_files(from_path=file_path, to_path='C:/cygwin64/home/wen_x', file_name='/SARS-CoV-2_nucleotide_fasta.txt')


# In[45]:


# move results.xml file back to my file_path for further analysis
move_files(from_path='C:/cygwin64/home/wen_x', to_path=file_path, file_name='/results.txt')


# In[44]:


# Q6. analyze blast result
analyze_blast_result2(file_path + '/results.txt')


# In[5]:


# additional analysis - tblastx result
analyze_blast_result2(file_path + '/tblastx_results.txt')


# In[ ]:




