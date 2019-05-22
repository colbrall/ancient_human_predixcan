'''
given an Ensembl or rs id, searches for PrediXcan models in a SQLite database

prints all models for the gene, or all models using the SNP, and what tissues
they're trained in

USAGE
python predixcan_search.py ID PATH/TO/DATABASE
ID can be an rsID or an Ensembl ID
use * to search multiple databases
'''

import commands
import sys
import string
import sqlite3

def geneSet(files,snps):
    f = [line.strip() for line in open(snps, 'r')]
    gene_set = set([])
    for item in files:
        if not item.endswith(".db"): continue
        query = "sqlite3 " + item
        for snp in f:
            s = str.split(snp,"\t")[7]
            sql_q = query + ' "SELECT gene FROM weights WHERE rsid=%s' % ("'" + s + "'" + '"')
            output = commands.getstatusoutput(sql_q)[1]
            if len(output) > 0:
                for gene in str.split(output,"\n"):
                    gene_set.add(gene)
    for gene in gene_set:
        print(gene)

def geneDist(files,snps):
    f = [line.strip() for line in open(snps, 'r')]
    print("File\tGene\tN_missing\tTotal")
    for item in files:
        gene_dict = {}
        if not item.endswith(".db"): continue
        cursor = sqlite3.connect(item).cursor()
        for snp in f:
            s = snp.strip()
            sql_snp = 'SELECT gene FROM weights WHERE rsid=%s' % ("'" + s + "'")
            for gene in cursor.execute(sql_snp):
                if gene[0] in gene_dict.keys():
                    gene_dict[gene[0]] = (gene_dict[gene[0]][0]+1, gene_dict[gene[0]][1])
                else:
                    sql_q = 'SELECT "n.snps.in.model" FROM extra WHERE gene=%s' % ('"' + gene[0] + '"')
                    output = cursor.execute(sql_q).fetchone()[0]
                    gene_dict[gene[0]] = (1,output)
        cursor.close()
        for gene in gene_dict.keys():
            print("%s\t%s\t%s\t%s" % (item,gene,gene_dict[gene][0],gene_dict[gene][1]))

def main():
    s = sys.argv[1]
    paths = sys.argv[2:]
    if s.startswith('rs'):
        print '\nsearching for SNP ' + s + '....\n'
        sql_q = ' "SELECT * FROM weights WHERE rsid=%s' % ("'" + s + "'" + '"')
    elif s.startswith('ENS'):
        print '\nsearching for gene ' + s + '....\n'
        sql_q = ' "SELECT * FROM weights' + '"' + " | grep " + s
    else:
        geneSet(paths,s)
        #geneDist(paths,s)
        sys.exit()

    for path in paths:
        p = path.split('/')[-1].split('.')[0]
        query = "sqlite3 " + path

        output = commands.getstatusoutput(query+' "'+"PRAGMA table_info('weights')"+'"')[1]
        head = ""
        for line in output.split("\n"):
            head += line.split('|')[1] + "|"
        output = commands.getstatusoutput(query + sql_q)[1]
        if len(output) > 0:
            print "\n" + p # print current database name
            print head[:-1] #print header
            print output + "\n"
        else:
            print "NONE- " + p


if __name__ == "__main__":
    main()
