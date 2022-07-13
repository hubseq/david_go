#
# Modified to take in a CSV or TXT file of input gene IDs. Must be in a column named gene_id and MUST be Ensembl IDs for now.
#

#import ssl
#ssl._create_default_https_context = ssl._create_unverified_context

import sys, os
sys.path.append('../')

import logging
import traceback as tb
import suds.metrics as metrics
from tests import *
from suds import *
from suds.client import Client
from datetime import datetime
import pandas as pd
sys.path.append('/global_utils/src/')
import module_utils

def isType( v ):
        """ Checks type of variable v. Eventually move this to utils()
        """
        try:
                tmp = int(v)
                return 'int'
        except:
                try:
                        tmp = float(v)
                        return 'float'
                except:
                        return 'str'

                
def chartReport( arg_list ):
        """ Creates signifcant GO terms output given list of gene_ids. 
        Can also pass conditions for which gene_ids to use, based on values in other columns.
        -i input_files
        -o output_dir
        -cond conditionals - e.g., pvalue<0.05,log2FoldChange>1
              valid conditionals: <, >, <=, >=, ==, !=
        -name analysis_name
        """        
        print('ARG LIST: {}'.format(str(arg_list)))
        input_args = module_utils.getArgument( arg_list, '-i', 'list' )
        output_dir = module_utils.getArgument( arg_list, '-o', 'list' )
        go_analysis_name = module_utils.getArgument( arg_list, '-name' )
        conditionals = module_utils.getArgument( arg_list, '-cond', 'implicit', '' )        
        if output_dir != []:
	        os.chdir( output_dir[0] )
        if input_args == []:
                print('INPUT FILE NOT SPECIFIED')
                return -1
        if input_args[0].lower().endswith('.csv'):
                df = pd.read_csv( input_args[0] )
        else:
                df = pd.read_csv( input_args[0], sep='\t')
        df.dropna(inplace=True)

        # reduce dataframe based on conditionals if provided
        operators = ['<','>','=','==','<=','>=','!=']
        if conditionals not in ['', []]:
                for c in conditionals.split(','):
                        c = c.lstrip(' ').rstrip(' ')
                        if '<=' in c:
                                [col, val] = c.split('<=')
                                if isType(val) == 'float':
                                        df = df[df[col]<=float(val)]
                                elif isType(val) == 'int':
                                        df = df[df[col]<=int(val)]
                        elif '>=' in c:
                                [col, val] = c.split('>=')
                                if isType(val) == 'float':
                                        df = df[df[col]>=float(val)]
                                elif isType(val) == 'int':
                                        df = df[df[col]>=int(val)]
                        elif '<' in c:
                                [col, val] = c.split('<')
                                if isType(val) == 'float':
                                        df = df[df[col]<float(val)]
                                elif isType(val) == 'int':
                                        df = df[df[col]<int(val)]
                        elif '>' in c:
                                [col, val] = c.split('>')
                                if isType(val) == 'float':
                                        df = df[df[col]>float(val)]
                                elif isType(val) == 'int':
                                        df = df[df[col]>int(val)]
                        elif '==' in c:
                                [col, val] = c.split('==')
                                if isType(val) == 'float':
                                        df = df[df[col]==float(val)]
                                elif isType(val) == 'int':
                                        df = df[df[col]==int(val)]                                        
                                else:
                                        df = df[df[col]==val]
                        elif '!=' in c:
                                [col, val] = c.split('!=')
                                if isType(val) == 'float':
                                        df = df[df[col]!=float(val)]
                                elif isType(val) == 'int':
                                        df = df[df[col]!=int(val)]                                        
                                else:
                                        df = df[df[col]!=val]                                        
                        
                                        
        print('GENE IDS: {}'.format(str(list(df['gene_id']))))
        inputIds = ','.join(list(df['gene_id'])) if 'gene_id' in list(df.columns) else ''
        
        outFile = os.path.join( output_dir[0], 'davidgo.{}.txt'.format('goterms' if go_analysis_name in ['', []] else go_analysis_name) )
        
        errors = 0
        
        setup_logging()
        logging.getLogger('suds.client').setLevel(logging.DEBUG)

        url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    
        print('url={}'.format(url))

        #
        # create a service client using the wsdl.
        #
        client = Client(url)
        client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
        #
        # print the service (introspection)
        #
        #print client
        
        #authenticate user email 
        client.service.authenticate('jerry@hubseq.com')

        #add a list 
        # inputIds = 'ENSMUSG00000084989,ENSMUSG00000026418,ENSMUSG00000100627,ENSMUSG00000104302,ENSMUSG00000055872,ENSMUSG00000044921,ENSMUSG00000020415,ENSMUSG00000090173,ENSMUSG00000020553,ENSMUSG00000096349,ENSMUSG00000095590,ENSMUSG00000021098,ENSMUSG00000097276,ENSMUSG00000015533,ENSMUSG00000022123,ENSMUSG00000063234,ENSMUSG00000022805,ENSMUSG00000055811,ENSMUSG00000024421,ENSMUSG00000024552,ENSMUSG00000042064,ENSMUSG00000087588,ENSMUSG00000058174,ENSMUSG00000043621,ENSMUSG00000085395,ENSMUSG00000093975,ENSMUSG00000078496,ENSMUSG00000081301,ENSMUSG00000096440,ENSMUSG00000078160,ENSMUSG00000090326,ENSMUSG00000004415,ENSMUSG00000067714,ENSMUSG00000054966,ENSMUSG00000103887,ENSMUSG00000053111,ENSMUSG00000031070,ENSMUSG00000054320,ENSMUSG00000032494,ENSMUSG00000031170'

        # idType = 'AFFYMETRIX_3PRIME_IVT_ID'
        idType = 'ENSEMBL_GENE_ID'
        listName = 'make_up'
        listType = 0
        print(client.service.addList(inputIds, idType, listName, listType))

        #print client.service.getDefaultCategoryNames()
        # setCategories
        #categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))

        #getChartReport
        thd = 0.1
        ct = 2
        chartReport = client.service.getChartReport(thd,ct)
        chartRow = len(chartReport)
        print('Total chart records:{}'.format(str(chartRow)))

        #parse and print chartReport
        resF = outFile
        with open(resF, 'w') as fOut:
                fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
                for simpleChartRecord in chartReport:
                        categoryName = simpleChartRecord.categoryName
                        termName = simpleChartRecord.termName
                        listHits = simpleChartRecord.listHits
                        percent = simpleChartRecord.percent
                        ease = simpleChartRecord.ease
                        Genes = simpleChartRecord.geneIds
                        listTotals = simpleChartRecord.listTotals
                        popHits = simpleChartRecord.popHits
                        popTotals = simpleChartRecord.popTotals
                        foldEnrichment = simpleChartRecord.foldEnrichment
                        bonferroni = simpleChartRecord.bonferroni
                        benjamini = simpleChartRecord.benjamini
                        FDR = simpleChartRecord.afdr
                        rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(FDR)]
                        fOut.write('\t'.join(rowList)+'\n')
                print('write file: ', resF, ' finished!')
