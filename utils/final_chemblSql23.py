#!/usr/bin/python

# Usage: chemblsql.py chembl_target_file
#
# This script reads a file of lines of tab delimited:
# contigname   chembl_targetid
#
# Outputs fields from the chembl mysql database
# Tab delimited with headers
#
# Ross Hall and Neil Young, 2010
# Modified by PKo, 2014
# And again by Ross hall 26/11/2014 for Chembl19
# Modified 09.03.2018 by Andreas Stroehlein for ChEMBL 23

import os, sys, MySQLdb

argc = len(sys.argv)
if argc != 2: 
    print 'Usage: chemblsql.py chembl_target_file'
    exit(1);
	
parameters = sys.argv
infile = parameters[1]

try:
    conn = MySQLdb.connect (host = "localhost",user ="bio",passwd = "password",db = "chembl_23")
except MySQLdb.Error, e:
    print "Error %d: %s" % (e.args[0], e.args[1])
    sys.exit (1)

handle = None
try: 
    handle = open(infile, 'r')
except:
    print "Trouble opening " + infile
    sys.exit(-1)	
		
tarray = []
for line in handle:
	words = line.rstrip()
	tarray.append(words)
handle.close()

headlist = ["#Accession number",
			"ChEMBL Target ID",
            "Target type",
            "Number of Lipinski's Rule-of-5 violations",
            "Rule-of-three pass",
            "Target name",
            "Compound name",
            "Compound synonyms",
            "ChEMBL compound ID",
            "Canonical SMILES",
            "Drug phase",
            "Natural product?",
            "First approval",
            "Patent No.",
            "Patent expiry date",
            "Assay type",
            "Activity type",
            "Activity relation",
            "Activity unit",
            "Activity value"
            ]
				
print '\t'.join(headlist)		

c = 0
for acc_id in tarray:
    #print acc_id
    cursor = conn.cursor()
    cursor.execute("""SELECT DISTINCT
  component_sequences.accession,
  target_dictionary.chembl_id,
  target_dictionary.target_type,
  compound_properties.num_lipinski_ro5_violations,
  compound_properties.ro3_pass,
  target_dictionary.pref_name,
  compound_records.compound_name,
  molecule_synonyms.synonyms,
  molecule_dictionary.chembl_id,
  compound_structures.canonical_smiles,
  molecule_dictionary.max_phase,
  molecule_dictionary.natural_product,
  molecule_dictionary.first_approval,
  product_patents.patent_no,
  product_patents.patent_expire_date,
  assays.assay_type,
  activities.standard_type,
  activities.standard_relation,
  activities.standard_units,
  activities.standard_value

  FROM component_sequences
  JOIN target_components ON (target_components.component_id = component_sequences.component_id)
  JOIN target_dictionary ON (target_dictionary.tid = target_components.tid)
  JOIN assays ON (assays.tid = target_dictionary.tid)
  JOIN activities ON (activities.assay_id = assays.assay_id)
  JOIN compound_records ON (compound_records.record_id = activities.record_id)
  LEFT JOIN compound_properties ON (compound_properties.molregno = compound_records.molregno)
  JOIN molecule_dictionary ON (molecule_dictionary.molregno = compound_records.molregno)
  JOIN compound_structures ON (compound_structures.molregno = molecule_dictionary.molregno)
  LEFT JOIN molecule_synonyms ON (molecule_synonyms.molregno = molecule_dictionary.molregno)
  LEFT JOIN formulations ON (formulations.molregno = compound_records.molregno) # the FK here used to be record_id which doesn't work. molregno works fine
  LEFT JOIN products ON (products.product_id = formulations.product_id)
  LEFT JOIN product_patents ON (product_patents.product_id = formulations.product_id)
  WHERE component_sequences.accession=%s""", (acc_id,))

    while True:
        row = cursor.fetchone()
        #print row
        if row == None:
          #print row
          break
        
        for x in row:
          print str(x) + "\t",
        print
    cursor.close ()
    c += 1
conn.close ()