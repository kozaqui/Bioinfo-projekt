__author__ = 'Kozel'
# go through PDB database for the first time
import re
import sqlite3
from ftplib import FTP
import gzip
from datetime import datetime
import time
import csv

from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import ModifiedRNAData


start_time = time.time()

### create/connect to a SQLite database
connection = sqlite3.connect('modRNA.db')
c = connection.cursor()

'''
########## for testing, delete tables in db
c.execute("DROP TABLE PDB_checklist")
c.execute("DROP TABLE mod_RNA_entities")
connection.commit()
########################
'''

def checkTableExist(cursor,tableName):
    """
    Checks if a table with a given name exists in the db, and returns
    a Boolean
    """
    cursor.execute("SELECT COUNT(*) FROM sqlite_master "
                   "WHERE type='table' "
                   "AND name=?", [(tableName)]) # TODO check sqlite syntax
    if cursor.fetchone()[0] == 1:
        return True
    else:
        return False



if not checkTableExist(c, "PDB_checklist"):
    c.execute("CREATE TABLE "
              "PDB_checklist(PDB_id TEXT, mod_date DATE, no_of_modRNA INT)")

if not checkTableExist(c, "mod_RNA_entities"):
   c.execute("""CREATE TABLE
              mod_RNA_entities(base_id INTEGER PRIMARY KEY,
                               PDB_id TEXT, mono_abbreviation TEXT,
                               mono_name TEXT, strand_ID TEXT,
                               residue INT)""")

if not checkTableExist(c, "mod_RNA_list"):
   c.execute("""CREATE TABLE
             mod_RNA_list(abbreviation TEXT, full_name TEXT)""")
   # '''
   #  c.execute("""CREATE TABLE
   #            mod_RNA_entities(base_id INTEGER PRIMIARY KEY,
   #                             PDB_id TEXT, mol_name TEXT,
   #                             base_name TEXT, chain_ID TEXT,
   #                             residue INT, PDB_link TEXT,
   #                             file_path TEXT)""")
   #                             '''

# TODO create a table with modified RNA list from Modomics
'''
if not checkTableExist(c, "modRNA_to_PDB"):
    c.execute("CREATE TABLE "
              "modRNA_to_PDB(PDB_id TEXT, mod_date DATE, no_of_modRNA INT)")
'''
3
# connect to Protein Data Bank through FTP
ftp = FTP('ftp.wwpdb.org')
ftp.login(passwd="") # login=anonymous, by default
 # go to a directory with some modRNA
ftp.cwd('/pub/pdb/data/structures/divided/mmCIF/lu')
# ftp.cwd('/pub/pdb/data/structures/all/mmCIF')


filenames = ftp.nlst()

sql_date_format = "%Y-%m-%d"
ftp_date_format = "%Y%m%d"
canonical_bases = ['A','U','G','C']

no_to_test = 0


for filename in filenames:
    PDB_ID = filename[0:4] # using the fact that PDB ID is 4-symbols
    print(PDB_ID)
    mod_base_count = 0

    # check if already in checked db
    db_date = c.execute("SELECT mod_date FROM PDB_checklist "
                        "WHERE PDB_id = ?", [(PDB_ID)])

    # i.e. if there's already a record (mod_date list is not empty)
    if not db_date:
        # check the modification date on ftp server
        db_date = datetime.strptime(db_date[0][0], sql_date_format)
        ftp_date = datetime.strptime(ftp.sendcmd('MDTM ' + filename)[4:12],
                                     ftp_date_format)
        if ftp_date > db_date:
            # TODO Behave, when there is a record, but has a new modified date
            a = 2

    # there's no record yet
    else:
        # if not
        # ftp_date = datetime.strptime(ftp.sendcmd('MDTM ?', [(filename)])[4:],
        #                              ftp_date_format)

        ftp_date = ftp.sendcmd('MDTM ' + filename)
        ftp_date = datetime.strptime(ftp_date[4:12], ftp_date_format)

        fileGz = open('temporary.gz','wb')
        ftp.retrbinary('RETR ' + filename, fileGz.write)
        fileGz.close()

        fileCif = open('temporary.cif','wb')
        f = gzip.open('temporary.gz', 'rb')
        fileCif.write(f.read())
        f.close()
        fileCif.close()

        # TODO Extract the header from the mmCIF file (no info about coordinates)

        mmCIF = MMCIF2Dict('temporary.cif')
        print(mmCIF['_entity_poly.type'])

##################
        if isinstance(mmCIF['_entity_poly.type'], list):
            for i, polymer in enumerate(mmCIF['_entity_poly.type']):
                if polymer == 'polyribonucleotide':
                    poly_id = mmCIF['_entity_poly.entity_id'][i]

                    for j, id in enumerate(mmCIF['_entity_poly_seq.entity_id']):
                        # check if there are modified bases in given strans
                        if id == poly_id and mmCIF['_entity_poly_seq.mon_id'][j] not in canonical_bases:
                            # mmCIF['_entity_poly_seq.mon_id'][j] # TODO: go through list of monomers just once
                            #mmCIF['_entity_poly_seq.mod_id']= 8
                            # else no. = 0
                            # if not, no. = 0
                            # if updated proceed
                            mod_base_count += 1
                            mono_name = mmCIF['_entity_poly_seq.mon_id'][j]

                            c.execute("INSERT INTO mod_RNA_entities VALUES(1,?,?,?,?,?)",
                                      (PDB_ID, mono_name,
                                       mmCIF['_chem_comp.name'][mmCIF['_chem_comp.id'].index(mono_name)],
                                       mmCIF['_entity_poly.pdbx_strand_id'][i],
                                       mmCIF['_entity_poly_seq.num'][j]))
                            connection.commit()
        else:
            if mmCIF['_entity_poly.type'] == 'polyribonucleotide':
                poly_id = mmCIF['_entity_poly.entity_id']

                for j, id in enumerate(mmCIF['_entity_poly_seq.entity_id']):
                    # check if there are modified bases in given strans
                    if id == poly_id and mmCIF['_entity_poly_seq.mon_id'][j] not in canonical_bases:
                        # mmCIF['_entity_poly_seq.mon_id'][j] # TODO: go through list of monomers just once
                        #mmCIF['_entity_poly_seq.mod_id']= 8
                        # else no. = 0
                        # if not, no. = 0
                        # if updated proceed
                        mod_base_count += 1
                        mono_name = mmCIF['_entity_poly_seq.mon_id'][j]

                        c.execute("INSERT INTO mod_RNA_entities VALUES(NULL,?,?,?,?,?)",
                                  (PDB_ID, mono_name,
                                   mmCIF['_chem_comp.name'][mmCIF['_chem_comp.id'].index(mono_name)],
                                   mmCIF['_entity_poly.pdbx_strand_id'],
                                   mmCIF['_entity_poly_seq.num'][j]))
                        connection.commit()
    '''
    c.execute("""CREATE TABLE
              mod_RNA_entities(base_id INTEGER PRIMIARY KEY,
                               PDB_id TEXT, mono_abbreviation TEXT,
                               mono_name TEXT, strand_ID TEXT,
                               residue INT""")
    '''
    c.execute("INSERT INTO PDB_checklist VALUES(?,?,?)",
              (PDB_ID, ftp_date.strftime(sql_date_format), mod_base_count))
    connection.commit()



                    # '''
                    #     if 'polyribonucleotide' in mmCIF['_entity_poly.type']: #??
                    #     # if there is a polyribonucleic acid
                    #         # check if there's more than one strand
                    #         if mmCIF['_entity_poly_type'] > 1:
                    #
                    #             polyrna_indices = [i for i, x in enumerate(mmCIF['_entity_poly.type'])
                    #                                if x == "polyribonucleotide"]
                    #                                '''


    #                 else:
    #                 c.execute("INSERT INTO PDB_checklist VALUES(?,?,?)",
    #           (PDB_ID, ftp_date.strftime(sql_date_format), 0))
    # connection.commit()

    no_to_test += 1
    if no_to_test == 2:
    #    continue
        break

# leave FTP nicely
ftp.quit()

## read data from modomics db file
# csv


# TODO retrieve all not-duplicated PDB_ID records containing given modRNA

modRNA_data = c.execute("""SELECT DISTINCT mono_abbreviation, mono_name
                           FROM mod_RNA_entities""")

for abbrev, name in modRNA_data.fetchall():
    c.execute("INSERT INTO mod_RNA_list VALUES(?,?)",(abbrev,name))
connection.commit()

modRNA2PDBDict = {}

#modRNA_data = c.execute("""SELECT abbreviation FROM mod_RNA_list""")
#??????????????????????????????????????????????????????????????
modRNA_data = c.execute("""SELECT DISTINCT mono_abbreviation
                           FROM mod_RNA_entities""")
for abbrev in modRNA_data.fetchall():
    PDB_list = []
    c.execute("""SELECT DISTINCT PDB_id FROM mod_RNA_entities
                 WHERE mono_abbreviation = ?""", [(abbrev[0])])
    for record in c.fetchall():
        PDB_list.append(record[0])
    modRNA2PDBDict[abbrev] = PDB_list

print("Modified RNA to PDB structures dictionary:")
print(modRNA2PDBDict)

end_time = time.time()
elapsed_time = end_time - start_time
print("Time to proceed: ", elapsed_time)

# check for changes in PDB database

# create dictionary mapping modified RNAs to PDB structures