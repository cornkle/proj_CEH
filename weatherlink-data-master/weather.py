#!/usr/bin/env python

import optparse
import datetime
import httplib
import csv
import os
import sys
import struct
import itertools

import progdata


def optionSetup():
    """ 
    Use optparse module to setup the program's command line options
    and defaults. The funtion returns the tuple (options, args) 
    where options.filename is the provided filename, etc.
    """
    p = optparse.OptionParser()
    p.add_option("-f", "--file", dest="filename", 
                 default="%s/weatherlink.csv" % os.getenv("HOME"),
                 help="csv file to store data", metavar="FILE")
    p.add_option("-q", "--quiet",
                 action="store_false", dest="verbose", default=True,
                 help="don't print status messages to stdout")
    return p.parse_args()


def readSettings(filename):
    """
    Read the client's weatherlink.com account username, password,
    and list of fields to exclude. 
    """
    try:
        f = open(filename, 'r') # open in append binary mode
    except IOError as err:
        print(sys.stderr, "Issue opening file: %s", filename)
        print(sys.stderr, err)
        exit(1)
    with f:
        username = f.readline().strip()
        password = f.readline().strip()
        exclusions = f.read()
        exclusions = "".join(exclusions.split()).split(',') 
        fields = [f for f in progdata.ordFieldsA if f not in exclusions]
        return username, password, fields


def timestamp(t):
    """
    Takes a datetime.datetime object and returns corresponding integer
    necessary in request to retrieve data. This value is effectively
    just the two byte date stamp followed by the two byte time stamp.
    >>> import datetime
    >>> tstamp(datetime.datetime(2012, 6, 17, 15, 5)) # 6/17/12 15:05
    416351713
    >>> tstamp(None)
    0
    """
    if t == None:
        return 0
    vantageTimeStamp = (100 * t.hour) + t.minute
    vantageDateStamp = t.day + (t.month << 5) + ((t.year - 2000) << 9)
    return (vantageDateStamp << 16) + vantageTimeStamp


def initCSV(reader, writer, fields):
    """
    make sure that the CSV file has the appropriate headers -- if they
    are missing, then write them, if they do not match up with 
    the user's specified fields then output a warning message. Finally,
    return date and time of the most recent record as a 
    datetime.datetime object (or None if there are no records).
    """
    try:
        headers = reader.next()
    except StopIteration: # handle empty file
        print("Writing headers to csv.")
        writer.writerow(fields)
        return None
    else: 
        if headers != fields:
            if len(set(fields) - set(headers)): 
                print(sys.stderr, "Missing fields in csv header: %s" \
                    % list(set(fields) - set(headers)))
            if len(set(headers) - set(fields)):
                print(sys.stderr, "Excluded fields in csv header: %s" \
                    % list(set(headers) - set(fields)))
            if raw_input("Enter 'q' to abort.\n") == 'q': 
                exit(1)
    finally:
        oldest = None
        for row in reader:
            try:
                t = ' '.join(row[0:2])    
                t = datetime.datetime.strptime(t, "%m/%d/%y %H:%M")
            except:
                # Unsuccessful attempt to create datetime.datetime from the row
                # (could be an empty line, bogus data, etc.); skip to next row
                continue 
            else: 
                oldest = t
        return oldest


def retrieveData(conn, username, password, t):
    """
    returns an httplib.HTTPReponse object containing the requested
    data from weatherlink.com after time t. 
    """
    url = "/webdl.php?timestamp=" + str(t) + "&user=" + username + "&pass=" + \
          password + "&action="
    conn.request("GET", url + "headers") # request with action=headers
    res = conn.getresponse()

    if res.status != 200:
        print(sys.stderr, "Could not access weatherlink server.")
        print(sys.stderr, "HTTP status code: %s: %s" % res.status, res.reason)
        exit(1)
    # parse out number of records
    numRecords = res.read().split("\n")[1].split("=")[1]
    print("Preparing to import %s records..." % numRecords)

    conn.request("GET", url + "data")    # request with action=data
    res = conn.getresponse()
    if res.status != 200:
        print(sys.stderr, "Could not retrieve data from server.")
        print(sys.stderr, "HTTP status code: %s: %s" % res.status, res.reason)
        exit(1)
    return res, int(numRecords)


def appendRecord(writer, record, fields, mask):
    """
    append the appropriate data from the 52 byte record onto the 
    provided local CSV file. This function handles unpacking the
    data and interpreting the bytes.
    """
    if record == progdata.DASHED:
        print(sys.stderr, "Skipping dashed record.")
        return
    # unpack 52 byte record into 39 data values
    vals = struct.unpack(progdata.fmtA, record) 
    # delete vals listed in dotfile
    vals = list(itertools.compress(vals, mask)) 

    recordDict = {f : progdata.format_map[f](v) for f, v in zip(fields, vals)}
    writer.writerow(recordDict)


def retrieve_data(filename,username,password,fields):
    # Create mask to delete unnecessary data in appendRecord()  
    mask = [0 if field not in fields else 1 for field in progdata.ordFieldsA]
    
    try:
        f = open(filename, 'a+b') # open in append mode
    except IOError as err:
        print(sys.stderr, "Issue opening file: %s", filename)
        print(sys.stderr, err)
        exit(1)
    with f:
        # Handle CSV setup, return t (timestamp of most recent record on file)
        t = initCSV(csv.reader(f), csv.writer(f), fields) 
        # convert datetime.datetime timestamp to WeatherLink format
        t = timestamp(t)
    
        # Retrieve all records more recent than t
        conn             = httplib.HTTPConnection("weatherlink.com")
        res, numRecords  = retrieveData(conn, username, password, t)
    
        # Create writer to append CSV
        writer = csv.DictWriter(f, fields) 

        # Read each retrieved 52 byte record and append it to the CSV file
        for i in range(numRecords):
            record = res.read(52) 
            appendRecord(writer, record, fields, mask) 

        conn.close()

def main():
    (options, arguments)         = optionSetup()
    (username, password, fields) = readSettings("%s/.weatherlink" % os.getenv("HOME")) 
    retrieve_data(options.filename,username,password,fields)


if __name__ == '__main__':
    main()
