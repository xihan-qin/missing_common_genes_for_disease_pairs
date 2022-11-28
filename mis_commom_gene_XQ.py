

"""
    After running the script You will get a plot showing the roc score frequency
    and you can find intermediate results in the output folder "Data".

    To run the script, type the line below in the command prompt.
    python3 mis_commom_gene_XQ.py

    You can also modify the code, comment out certain steps.

"""

################################################################################
import sys
import os
import glob
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from collections import Counter
import networkx as nx
import numpy as np
import random
# functions from other paper
from separation import *
"""
functions include:
read_network, read_gene_list, remove_self_links,get_pathlengths_for_single_set,
get_pathlengths_for_two_sets, calc_single_set_distance, calc_set_pair_distances
"""

################################################################################
################################################################################
"""
    main body of the program
"""

def main():
    """
        you can comment out selective step if you already have some steps done
    """

    #------------------------Pre_setting_block_start---------------------------#
    # input folders
    input_folder = "input_files"

    # !!! Make sure you have the input files to work with
    network_file = f'{input_folder}/DataS1_interactome.tsv' # the whole PPI graph dataset
    Disease_gene_list_file = f'{input_folder}/DataS2_disease_genes.tsv' # the file contains gene lists for diseases
    Disease_pair_file = f'{input_folder}/DataS4_disease_pairs.tsv'  # the file stores the true disease pairs with their common genes and s_AB

    # Provide the paths and file name for intermediate results
    output_folder = 'Data'
    diseasepath = f'{output_folder}/Disease_genes' # stores the disease gene sets
    rocpath = f'{output_folder}/ROC'    # stores gene, s_AB_delta, label for each disease pair
    r = f'{output_folder}/R'    # stores the data to calculate roc for each disease pair
    negsab_file_path = f'{output_folder}/negsab.tsv'   # stores the disease pairs that have negative S_AB
    roc_file_path = f'{output_folder}/roc.txt'  # stores disease_pair, roc_auc, uncommon_gene_cont, common_gene_cont, s_AB_delta

    # if these folders not exist, create the path
    folder_path_list = [diseasepath, rocpath, r]
    check_folders_existence(folder_path_list)

    #------------------------Pre_setting_block_end-----------------------------#

    #------------------------main_steps_start----------------------------------#
    # you can selectively comment out the steps based on your case and need.

    # step1 -- generate files for disease gene lists
    generate_disease_gene_list_files(Disease_gene_list_file, diseasepath)

    # step2 -- generate file 'negsab.tsv'
    generate_negsab(Disease_pair_file, negsab_file_path)

    # step3 -- select how many disease pairs to go over with
    no_disease_pairs = 5 # e.g. 5 to go over 5 pairs; 'all' to go over all pairs
    disease_pair_list, disease_pair_sab_dict = get_disease_pair_list(no_disease_pairs, negsab_file_path) # parse 'negsab.tsv'

    # step4 -- get the network graph
    G, all_genes_in_network = get_graph(network_file)

    # step5 -- Iteratively calculate s_AB_remove for common_gene and s_AB_add for uncommon_gene
    s_AB_for_all_genes(disease_pair_list, diseasepath, rocpath, G, all_genes_in_network)     # results are stored into files, to run all the paires takes several hours

    # step6 -- generate results for roc calculation -- s_AB_delta is calculated in this step
    roc_file(rocpath,r, disease_pair_sab_dict)       # results are stored into files

    # step7 -- calculate roc score
    roc_score(r, roc_file_path)    # results are stored into files

    # step8 -- plot histogram of roc_score vs frequency in disease_pairs
    hist_plot(roc_file_path)
    #------------------------main_steps_end------------------------------------#
################################################################################
#####--------------------------functions_start-----------------------------#####
def generate_disease_gene_list_files(Disease_gene_list_file, diseasepath):
    """
        generate files
        Each file contains gene list for a disease including GWAS and OMIM genes
        Ensure the generated files are the same by
        comparing with the provided gene list files from the Author
    """
    f = open(Disease_gene_list_file, 'r')
    for line in f:
        if line[0]=='#':
            continue
        row = line.strip().split('\t')
        disease_name = row[0]
        # if the length is 5, then there is just one gene list, either OMIM or GWAS gene list
        if len(row) == 5:
            disease_gene_list = row[4].split(";")
        # if the length is 6, then there are two gene lists: OMIM and GWAS,
        # there could be duplicate between two lists, need to remove duplicates
        if len(row) == 6:
            disease_gene_list = list(set(row[4].split(";") + row[5].split(";")))
        file_name_path = f'{diseasepath}/{disease_name}.txt'
        save_list_to_file(disease_gene_list,file_name_path)


#------------------------------------------------------------------------------#
def check_folders_existence(folder_path_list):
    """
        given a list of folder path, if not exist, create the path
    """
    for folder_path in folder_path_list:
        if not os.path.exists(folder_path):
        	os.makedirs(folder_path)
#------------------------------------------------------------------------------#
def save_list_to_file(given_list,file_name):
    """
        write each item into a file as a single line
    """
    f = open(file_name,"w")
    for item in given_list:
        f.write(f'{item}\n')
    f.close()

# mostlikely not needed, it needed, add these lines
# save_list_to_file(common_genes,"common_genes.txt")
# save_list_to_file(uncommon_genes,"uncommon_genes.txt")
#------------------------------------------------------------------------------#
def write_list_to_line_in_file(given_list,file):
    """
        write list into a line in a given file
        line:
        item_1 \t item_2 \t ... \t item__n \n
    """
    line = "\t".join([str(i) for i in given_list])
    file.write(f'{line}\n')

#------------------------------------------------------------------------------#
def generate_negsab(Disease_pair_file, negsab_path):
    """
        process "DataS4_disease_pairs"
        return a file contains disease pairs that have negative SAB only
    """
    f1 = open(Disease_pair_file, 'r')
    f2 = open(negsab_path, 'w')
    for line in f1:
        if line[0]=='#':
            continue
        row = line.strip().split('\t')
        (disease_a, disease_b, s_AB) = row[:3]
        if float(s_AB) < 0:
            write_list_to_line_in_file([disease_a, disease_b, s_AB],f2)
    f1.close()
    f2.close()

#------------------------------------------------------------------------------#
def get_disease_pair_list(no_disease_pairs, negsab_path):
    negsab_file = open(negsab_path, 'r')
    disease_pair_list = []
    disease_pair_sab_dict = {}
    pair_count = 0
    for line in negsab_file:
        pair_count += 1
        disease_a, disease_b, s_AB = line.strip().split("\t")
        disease_pair_list.append([disease_a, disease_b])
        disease_pair_sab_dict[f'{disease_a}&{disease_b}'] = s_AB
        if pair_count == no_disease_pairs:
            break
    return disease_pair_list, disease_pair_sab_dict

#------------------------------------------------------------------------------#
def get_graph(network_file):
    G  = read_network(network_file)
    # get all genes ad remove self links
    all_genes_in_network = set(G.nodes())
    remove_self_links(G)
    return G, all_genes_in_network

#------------------------------------------------------------------------------#
def s_AB_for_all_genes(disease_pair_list, diseasepath, rocpath, G, all_genes_in_network):
    for (disease_a,disease_b) in disease_pair_list:
        #separate common and uncommon genes from disease pair
        common_genes,uncommon_genes=find_genes(disease_a,disease_b,diseasepath)
        if len(common_genes) != 0 and len(uncommon_genes) != 0:
            s_AB_uncommon_genes(disease_a,disease_b,rocpath, uncommon_genes, diseasepath, G, all_genes_in_network)
            s_AB_common_genes(disease_a,disease_b,rocpath, common_genes, diseasepath, G, all_genes_in_network)

#------------------------------------------------------------------------------#
def find_genes (gene_file_1,gene_file_2,diseasepath):
    ''''subroutine to separate common and uncommon genes from disease pair'''
    f1=open(diseasepath+"/"+gene_file_1+".txt","r")
    f2=open(diseasepath+"/"+gene_file_2+".txt","r")
    gene_list_1 = []
    gene_list_2 = []
    for line in f1:
        gene = line.strip()
        if gene not in gene_list_1:
            gene_list_1.append(gene)
    for line in f2:
        gene = line.strip()
        if gene not in gene_list_2:
            gene_list_2.append(gene)
    common_genes= list(set([i for i in gene_list_1 if i in gene_list_2]))
    all_genes = set(gene_list_1+gene_list_2)
    uncommon_genes= [item for item in all_genes if item not in common_genes]
    return common_genes, uncommon_genes

#------------------------------------------------------------------------------#
def s_AB_calculation(G, gene_set_A, gene_set_B):
	# distances WITHIN the two gene sets:
    d_A = calc_single_set_distance(G,gene_set_A)
    d_B = calc_single_set_distance(G,gene_set_B)
	# distances BETWEEN the two gene sets:
    d_AB = calc_set_pair_distances(G,gene_set_A,gene_set_B)
	# calculate separation
    s_AB = d_AB - (d_A + d_B)/2
    return s_AB

#------------------------------------------------------------------------------#
def get_gene_set(diseasepath,disease,all_genes_in_network):
	# read gene set 1
    gene_file_path = diseasepath+"/"+disease+".txt"
    genes_set_full = read_gene_list(gene_file_path) # returns a set of genes in the file
	# removing genes that are not in the network:
    genes_set = genes_set_full & all_genes_in_network

    return genes_set

#------------------------------------------------------------------------------#
def s_AB_uncommon_genes (disease_a,disease_b,rocpath,uncommon_genes, diseasepath, G, all_genes_in_network):
    """
        calculate the s_ab if the uncommon gene is the shared the gene
        from the uncommon_genes, choose randomly max len 10 of the gene set a
        each gene in a is added to both disease_a set and disease_b set as a "shared" gene
        then calculate the s_ab score if this gene becomes the shared gene between these two sets
        and store s_ab scores in the files
    """
    fw=open(rocpath+"/"+disease_a+"&"+disease_b+".txt","w")
    a=[] # a list of uncommon genes, if more than 10, randomly choose 10 from the list
    if len(uncommon_genes) <= 10:
        a=uncommon_genes
    else:
        a = np.random.choice(uncommon_genes, size=10)

    genes_A = get_gene_set(diseasepath,disease_a, all_genes_in_network)
    genes_B = get_gene_set(diseasepath,disease_b, all_genes_in_network)
    com_id = 0  # common gene indicator. if common_gene, com_id =1; otherwise, com_id = 0

    for gene in a:
        # making uncommon gene a common gene for gene A, union the gene to genes_A
        genes_A_and_i =  genes_A | set(gene)
		# making uncommon gene a common gene for gene B
        genes_B_and_i  = genes_B | set(gene)

        s_AB_add_gene = s_AB_calculation(G, genes_A_and_i, genes_B_and_i)

		# save results:
        write_list_to_line_in_file([gene, s_AB_add_gene, com_id],fw)
    fw.close()

#------------------------------------------------------------------------------#
def s_AB_common_genes (disease_a,disease_b,rocpath,common_genes, diseasepath, G, all_genes_in_network):
    """
        from the uncommon_genes, choose randomly max len 5 of the gene set a
        original code:
        for each gene in a
            remove from disease disease_gene_set a, calculate the s_AB
            remove from disease disease_gene_set b, calculate the s_AB
            then, for each gene, there will be two lines,
            one for s_ab when removing from disease_gene_set_a,
            one for s_ab when removing from disease_gene_set_b
        This function is now changed to:
        for each gene in a
            remove gene from both sets and then calculate s_AB
    """
    fw=open(rocpath+"/"+disease_a+"&"+disease_b+".txt","a+")
    a=[] # a list of common genes, if more than 5, randomly choose 5 from the list
    if len(common_genes) <= 5:
        a=common_genes
    else:
        a = np.random.choice(common_genes, size=5)

    genes_A = get_gene_set(diseasepath,disease_a, all_genes_in_network)
    genes_B = get_gene_set(diseasepath,disease_b, all_genes_in_network)

    com_id = 1 # common gene indicator. if common_gene, com_id =1; otherwise, com_id = 0

    for gene in a:
        # the common_gene may not in the network thus removed
        # need to check if common_gene in both sets
        if (gene in genes_A) and (gene in genes_B):
            genes_A_no_i =  genes_A.copy()
            genes_A_no_i.remove(gene)

            genes_B_no_i =  genes_B.copy()
            genes_B_no_i.remove(gene)

            s_AB_rm_gene = s_AB_calculation(G, genes_A_no_i, genes_B_no_i)

    		# save results
            write_list_to_line_in_file([gene, s_AB_rm_gene, com_id],fw)
    fw.close()

#------------------------------------------------------------------------------#
def roc_file (rocpath, r, disease_pair_sab_dict):
    path = rocpath+'/*.txt'
    files=glob.glob(path)
    for file in files:
        file_name_alone = file.split('/')[-1].strip(".txt") # disease_a&disease_b e.g.abnormalities, multiple&chromosome disorders
        if file_name_alone in disease_pair_sab_dict:
            s_AB_orig = disease_pair_sab_dict[file_name_alone]
            f = open(file, 'r')
            fw= open(r+"/"+file_name_alone+".txt","w")
            for line in f:
                (gene, s_AB, com_id) = line.strip().split()
                # calculate s_AB_delta = s_AB_without_gene - s_AB_with_gene
                # if com_id is 1, then s_AB is the score when a comman gene is removed
                # if com_id is 0, then s_AB is the score when a comman gene is added
                if com_id == "1":
                    s_AB_delta = float(s_AB)-float(s_AB_orig)
                if com_id == "0":
                    s_AB_delta = float(s_AB_orig)-float(s_AB)
                #In original code: fw.write(gene +"\t"+ s_AB +"\t" + str(f1[2]) +"\t"+str(f1[3])+"\t"+str(float(x[1])-float(f1[3]))+"\t"+str(x[2])+"\n")
                #not sure what f1[2] is. f1[3] is assumed to be the S_ab for the disese pair in DataS4_disease_pairs file
                write_list_to_line_in_file([gene, s_AB, s_AB_orig, s_AB_delta, com_id],fw)
            f.close()
            fw.close()

#------------------------------------------------------------------------------#
def roc_score (r, roc_file_path):
    path = r+'/*.txt'
    files=glob.glob(path)
    fw=open(roc_file_path,"w")
    #fw.write("#Disease Pair"+"\t"+"ROC Score"+"\t"+"Zero Count"+"\t"+"One Count"+"\t"+"Comorbidity Value\n")
    for file in files:
        file_name_alone = file.split('/')[-1].strip(".txt") # e.g.abnormalities, multiple&chromosome disorders
        (dis_a, dis_b) = file_name_alone.split("&")

        f=open(file, 'r')
        data=[]
        res=[]
        for line in f:
            (gene, cal_s_AB, s_AB_orig, s_AB_delta, com_id)=line.strip().split("\t")
            data.append(float(s_AB_delta))
            res.append(int(com_id))
        f.close()
        ab= set(res)
        if len(ab) == 2:
            label_counter_dict=Counter(res) # {key:Counter} e.g.{0: 10, 1: 5} for [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
            false_positive_rate, true_positive_rate, thresholds = roc_curve(res, data)
            roc_auc = auc(false_positive_rate, true_positive_rate)
            ax=roc_auc_score(res, data)
            # the original code write cal_s_AB into file
            # I changed to s_AB_delta to write into the same position
            write_list_to_line_in_file([file_name_alone, roc_auc, label_counter_dict[0], label_counter_dict[1], s_AB_delta],fw)
    fw.close()

#------------------------------------------------------------------------------#
def hist_plot (roc_file_path):
    f = open(roc_file_path, 'r')
    x=[]
    for line in f:
        (disease_pair, roc_auc, uncommon_cont, common_cont, s_AB_delta)=line.strip().split("\t")
        x.append(float(roc_auc))
	#print x
    plt.hist(x)
    plt.title("ROC Score Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.show()
#####--------------------------functions_end-------------------------------#####
################################################################################
if __name__ == '__main__':
    main()

sys.exit(0)
