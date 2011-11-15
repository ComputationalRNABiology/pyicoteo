# coding=utf-8
import sys
import MySQLdb
from MySQLdb import OperationalError, ProgrammingError

def retrieve_chrinfo(schema):
        try:
            print 'trying to retrieve chromosome length from ucsc (%s)...'%schema
            db = MySQLdb.connection(host="genome-mysql.cse.ucsc.edu",user="genome",db=schema)
            db.query("""select chrom, size from chromInfo""")
            rows = db.store_result()
            outfile = file(schema, 'w')
            for row in rows.fetch_row(10000000):
                outfile.write('%s\t%s\n'%(row[0], row[1]))
            print 'done!'
        except (OperationalError, ProgrammingError):
            print 'The schema %s doesnt have the table chromInfo'%schema
        finally:
            db.close()

def main():
    schema = sys.argv[1]
    print schema
    #get all the database names from the information schema in ucsc
    if schema == 'all':
        print 'retrieving all schema names...'
        try:
            db = MySQLdb.connection(host="genome-mysql.cse.ucsc.edu",user="genome",db="information_schema")
            db.query("""select schema_name from schemata""")
            rows = db.store_result()
            schemas = []
            for row in rows.fetch_row(10000000):
                schemas.extend(row)

        finally:
            db.close()
        #retrieve all chrinfo that exist
        for schema in schemas:
            retrieve_chrinfo(schema)
    else:
        retrieve_chrinfo(schema)

main()
