__author__ = 'Kozel'
import sqlite3
from ftplib import FTP
import re
import os


def checkTableExist(cursor,tableName):
    """
    Checks if a table with a given name exists in the database indicated by
    the cursor, and returns a Boolean
    """
    cursor.execute("SELECT COUNT(*) FROM sqlite_master "
                   "WHERE type='table' "
                   "AND name=?", (tableName,))
    if cursor.fetchone()[0] == 1:
        return True
    else:
        return False

def populate_mod_RNA_table(sqlite_connection):
    # if components.cif is not present download it through ftp
    if not os.path.isfile("components.cif"):
        ftp = FTP('ftp.wwpdb.org')
        ftp.login(passwd="") # login=anonymous, by default
        ftp.cwd('/pub/pdb/data/monomers')
        ftp.retrbinary('RETR components.cif', open('components.cif').write)
        ftp.quit()


    #### connect to db
    # connection = sqlite3.connect('modRNA_test.db')
    c = sqlite_connection.cursor()

    table_existed = 1
    if not checkTableExist(c, "mod_RNA_list"):
       c.execute("""CREATE TABLE
                 mod_RNA_list(mono_id VARCHAR(3), full_name TEXT)""")
       table_existed = 0



    # regex for extracting chemical compund id and its type
    fields_of_interest_re = '_chem_comp\.id\s*(\S*).*?_chem_comp\.name\s*"?(.*?)"?\s*?\n.*?_chem_comp\.type.\s*"?(.*?)"?\s*?\n'

    monomers = re.findall(fields_of_interest_re, open("components.cif").read(), re.DOTALL)

    # if populating table for the first time
    if not table_existed:
        # record[0] = id, record[1] = full_name, record[2] = type
        for record in monomers:
            # there are possible values 'RNA LINKING' and 'L-RNA LINKING'
            if record[2].endswith('RNA LINKING'):
                # populate the table
                print("Record: ", record)
                c.execute('''INSERT INTO mod_RNA_list(mono_id, full_name)
                             VALUES(?,?)''',
                          (record[0],record[1]))
                # delete AGUC
                c.execute('''DELETE FROM mod_RNA_list
                             WHERE mono_id = 'A' OR mono_id = 'G'
                             OR mono_id = 'U' OR mono_id = 'C' ''')

    sqlite_connection.commit()