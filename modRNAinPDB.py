__author__ = 'Kozel'
# go through PDB database for the first time
import re
import sqlite3
from ftplib import FTP
import gzip
from datetime import datetime
import time
from io import StringIO
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
    Checks if a table with a given name exists in the database indicated by
    the cursor, and returns a Boolean
    """
    cursor.execute("SELECT COUNT(*) FROM sqlite_master "
                   "WHERE type='table' "
                   "AND name=?", (tableName,)) # in sqlite execute() syntax, argument is a tuple or list of tuples
    if cursor.fetchone()[0] == 1:
        return True
    else:
        return False


dictCheckedValue = lambda name, dictionary:\
    dictionary[name] if name in dictionary.keys() else "NA"
# this function checks if a given dictionary has a given key (name), if yes
# it returns the value for that key, it returns "NA" otherwise


def extract_mmCIF_header(ftp, filename):
    """
    this function extracts mmCIF header - i.e. information before atom
    coordinates - from a file downloaded via ftp in compressed .cif.gz format,
    and returns it as a python dictionary
    :param ftp:
    :param filename:
    :return:
    """
    fileGz = open('temporary.gz','wb')
    ftp.retrbinary('RETR ' + filename, fileGz.write)
    fileGz.close()

    fileCif = open('temporary.cif','wb')
    f = gzip.open('temporary.gz', 'rb')
    fileCif.write(f.read()) # write unzipped file
    f.close()
    fileCif.close()

    # Extract the header from the mmCIF file (no info about coordinates)
    fileCif = open('temporary.cif', 'r')
    fileHeader = open('header.cif','w')
    full_cif = fileCif.read()
    fileCif.close()
    # no relevant info after this field: _struct_biol.id
    fileHeader.write(full_cif[:full_cif.find("_struct_biol.id")])
    fileHeader.close()

    mmCIF = MMCIF2Dict('header.cif')
    return mmCIF


def insert_modRNA_records(c, mmCIF):
    """
    Check if there are modified RNAs, if tes update mod_RNA_entities table
    and update PDB_checklist table, if it's the first time
    :param c: SQLite cursor to the database
    :param mmCIF: python dictionary with structure info
    """
    mod_base_count = 0
    if '_entity_poly.type' in mmCIF.keys():
        print(mmCIF['_entity_poly.type']) # TODO delete this line after testing

    ##################
        # if there is more than one polymer
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

                            # insert modified RNA record to database
                            c.execute("INSERT INTO mod_RNA_entities VALUES(NULL,?,?,?,?,0)",
                                      (PDB_ID, mono_name,
                                       mmCIF['_entity_poly.pdbx_strand_id'][i],
                                       mmCIF['_entity_poly_seq.num'][j]))
        #                    connection.commit()
        # if there's only one polymer
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

                        # insert modified RNA record to database
                        c.execute("INSERT INTO mod_RNA_entities VALUES(NULL,?,?,?,?,0)",
                                  (PDB_ID, mono_name,
                                   mmCIF['_entity_poly.pdbx_strand_id'],
                                   mmCIF['_entity_poly_seq.num'][j]))
            #            connection.commit()
        '''
        c.execute("INSERT INTO PDB_checklist VALUES(?,?,?,?,?)",
                  (PDB_ID, ftp_date.strftime(sql_date_format),
                   mmCIF["_struct_keywords.pdbx_keywords"],
                   str(dictCheckedValue("_entity_src_nat.pdbx_organism_scientific", mmCIF)), # what if more than one organism? then use str()
                   mmCIF["_struct.title"]))
        # connection.commit()
        '''


# Create three tables in the database, if they don't exist yet

if not checkTableExist(c, "PDB_checklist"):
    c.execute("CREATE TABLE "
              "PDB_checklist(PDB_id TEXT, mod_date DATE, molecule_type TEXT,"
                            "organism_name TEXT, description TEXT)")

if not checkTableExist(c, "mod_RNA_entities"):
   c.execute("""CREATE TABLE
              mod_RNA_entities(base_id INTEGER PRIMARY KEY,
                               PDB_id TEXT, mono_abbreviation TEXT,
                               strand_ID TEXT, residue INT,
                               obsolete INT)""") # no boolean in sqlite

if not checkTableExist(c, "mod_RNA_list"):
   c.execute("""CREATE TABLE
             mod_RNA_list(abbreviation TEXT, full_name TEXT)""")

# TODO create a table with modified RNA list from Modomics

# TODO pub/pdb/data/status - check which structures have been addded/modified/obsolete
'''
if not checkTableExist(c, "modRNA_to_PDB"):
    c.execute("CREATE TABLE "
              "modRNA_to_PDB(PDB_id TEXT, mod_date DATE, no_of_modRNA INT)")
'''
# TODO ask user which part of db they want to search (all, subset(e.g. "aa")
# or a specific structure (e.g. 1ab3)
# TODO Write a whole new program? /pub/pdb/data/monomers/components.cif

# TODO function to go to a specified location on PDB FTP server
# connect to Protein Data Bank through FTP
ftp = FTP('ftp.wwpdb.org')
ftp.login(passwd="") # login=anonymous, by default
# go to a directory with some modRNA
ftp.cwd('/pub/pdb/data/structures/divided/mmCIF/lu')
# ftp.cwd('/pub/pdb/data/structures/all/mmCIF')




filenames = ftp.nlst()

sql_date_format = "%Y-%m-%d"    # e.g. 2015-01-27
ftp_date_format = "%Y%m%d"      # e.g. 20150127
canonical_bases = ['A','U','G','C']

no_to_test = 0

########################
## THE MAIN LOOP #######
########################


for filename in filenames:
    PDB_ID = filename[0:4] # using the fact that PDB ID is 4-symbols
    print(PDB_ID)

    # check if already in checked db
    c.execute("SELECT mod_date FROM PDB_checklist "
              "WHERE PDB_id = ?", (PDB_ID,))
    db_date = c.fetchone()
    print(db_date)
    # extract the modification date, from ftp server
    ftp_date = ftp.sendcmd('MDTM ' + filename)
    ftp_date = datetime.strptime(ftp_date[4:12], ftp_date_format)

    # i.e. if there's already a record (mod_date list is not empty)
    if db_date:
        # check the modification date on ftp server
        db_date = datetime.strptime(db_date[0], sql_date_format)
        if ftp_date > db_date:
            mmCIF = extract_mmCIF_header(ftp, filename)
            c.execute("""UPDATE mod_RNA_entities
                         SET obsolete = 1
                         WHERE PDB_id = ? AND obsolete = 0""", (PDB_ID, ))
            #insert new modRNA records
            insert_modRNA_records(c, mmCIF)
            c.execute("""UPDATE PDB_checklist
                         SET mod_date = ?
                         WHERE PDB_id = ?""",
                      (ftp_date.strftime(sql_date_format), PDB_ID))
            connection.commit()

    # there's no record yet
    else:
        # if not
        # ftp_date = datetime.strptime(ftp.sendcmd('MDTM ?', [(filename)])[4:],
        #                              ftp_date_format)

        mmCIF = extract_mmCIF_header(ftp, filename)
        insert_modRNA_records(c, mmCIF)
        c.execute("INSERT INTO PDB_checklist VALUES(?,?,?,?,?)",
                  (PDB_ID, ftp_date.strftime(sql_date_format),
                   mmCIF["_struct_keywords.pdbx_keywords"],
                   str(dictCheckedValue("_entity_src_nat.pdbx_organism_scientific", mmCIF)), # what if more than one organism? then use str()
                   mmCIF["_struct.title"]))
        connection.commit()

    no_to_test += 1
    if no_to_test == 15:
        continue
    #    break


# leave FTP nicely
ftp.quit()

## read data from modomics db file
# csv


# TODO retrieve all not-duplicated PDB_ID records containing given modRNA

mod_RNA_data = c.execute("""SELECT DISTINCT mono_abbreviation
                           FROM mod_RNA_entities""")

fileCompoundDict = open('Components-pub.cif','r')
mmCIF_compounds = fileCompoundDict.read()

for abbrev in mod_RNA_data.fetchall():
    c.execute("SELECT * FROM mod_RNA_list "
               "WHERE abbreviation = ?", (abbrev[0],))
    list_present = c.fetchone()

    if not list_present:
        start = mmCIF_compounds.find("data_"+abbrev[0])
        stop = mmCIF_compounds.find("_chem_comp.type",start)
        small_dict = MMCIF2Dict(StringIO(mmCIF_compounds[start:stop]))
        c.execute("INSERT INTO mod_RNA_list VALUES(?, ?)",(abbrev[0],small_dict['_chem_comp.name']))

connection.commit()


'''
    db_date = c.execute("SELECT mod_date FROM PDB_checklist "
                        "WHERE PDB_id = ?", (PDB_ID,))

    # extract the modification date, from ftp server
    ftp_date = ftp.sendcmd('MDTM ' + filename)
    ftp_date = datetime.strptime(ftp_date[4:12], ftp_date_format)

    # i.e. if there's already a record (mod_date list is not empty)
    if not db_date:
    '''

# Now create a python dictionary mapping modified RNAs to PDB structures
modRNA2PDBDict = {}
mod_RNA_data = c.execute("""SELECT abbreviation
                           FROM mod_RNA_list""")

for abbrev in mod_RNA_data.fetchall():
    PDB_list = []
    c.execute("""SELECT DISTINCT PDB_id FROM mod_RNA_entities
                 WHERE mono_abbreviation = ?""", (abbrev[0],))# before [(abbrev[0])])
    for record in c.fetchall():
        PDB_list.append(record[0])
    modRNA2PDBDict[abbrev] = PDB_list

print("Modified RNA to PDB structures dictionary:")
print(modRNA2PDBDict)

end_time = time.time()
elapsed_time = end_time - start_time
print()
print("Time the program was running:", elapsed_time ,"s.")