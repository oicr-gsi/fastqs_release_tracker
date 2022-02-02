# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 13:26:14 2022

@author: rjovelin
"""

import argparse
import subprocess
import os
from datetime import datetime
import time
from dateutil.relativedelta import relativedelta
import requests

   

def get_workflow_id(file_path):
    '''
    (str) -> int
    
    Parameters
    ----------
    
    Returns the workflow id from the file path or -1 if workflow id not in file path
           
    - file_path (str): File path
    '''
    
    workflow = -1
    
    k = file_path.split('/')
    for j in k:
        if j.isdigit():
            workflow = int(j)
    return workflow
    

def convert_to_epoch(date):
    '''
    (str) -> int
    
    Returns a date in epoch time
    
    Parameters
    ----------
    - date (str): Date formatted as Y.m.d H:M:S
    '''
    
    p = '%Y.%m.%d %H:%M:%S'
    date = int(time.mktime(time.strptime(date, p)))
    return date


def format_date(date):
    '''
    (str) -> str
    
    Returns date in the format Y.m.d H:M:S
    
    Parameters
    ----------
    - date (str): A date in the format 'Y-m-d H:M:S.MS'
    '''

    date = date[:date.index('.')].split()
    date = '.'.join(date[0].split('-')) + ' ' + date[1]
    return date


def collect_records_from_FPR(provenance, workflow):
    '''
    (str, str)
    
    Returns a list of records from File Provenance Report for a given workflow
    
    Paramaters
    ----------
    - provenance (str): Path to the file provenance report
    - workflow (str): Name of the workflow
    '''
    
    try:
        records = subprocess.check_output('zcat {0} | grep {1}'.format(provenance, workflow), shell=True).decode('utf-8').rstrip().split('\n')
    except:
        records = []
    finally:
        return records



def map_instrument_type(sequencer):
    '''
    (str) -> str

    Returns a generic intrument name for the sequencer model
    
    Parameters
    ----------
    - sequencer (str): Name of the instrument on which libraries are sequenced
    '''
    
    instrument = ''
    if 'miseq' in sequencer.lower():
        instrument = 'MiSeq' 
    elif 'nextseq' in sequencer.lower():
        instrument = 'NextSeq'
    elif 'hiseq' in sequencer.lower():
        instrument = 'HiSeq'
    elif 'novaseq' in sequencer.lower():
        instrument = 'NovaSeq'
    return instrument






def extract_fastqs(provenance, time_interval, keep_novaseq):
    '''
    (str, int, bool) -> (dict)
  
    Returns a dictionary with file path and file swid organied by project and run
            
    Parameters
    ----------
    - provenance (str): Path to File Provenance Report
    - time_interval (int): Number of months prior the current date from which records are considered
    - keep_novaseq (bool): Keep only novaseq runs if True
    '''
    
    # create a dict {project: {run: {filename: {'filepath': filepath, 'swid': swid}}}}
    D = {}
    
    # make a list of records from file provenance record
    records = collect_records_from_FPR(provenance, 'casava') + collect_records_from_FPR(provenance, 'bcl2fastq')
            
    # compute date from which records are considered and format it to Y.m.d H:M:S
    start_date = (datetime.now() - relativedelta(months=+abs(time_interval))).strftime('%Y.%m.%d %H:%M:%S')
    # convert to epoch time
    start_date = convert_to_epoch(start_date)
    
    # parse the records
    for i in records:
        i = i.rstrip().split('\t')
        # get file path and swid
        file_path, swid = i[46], int(i[44])
        # get filename
        filename = os.path.basename(file_path)
        # get project and run
        project, run_id = i[1], i[18]
        # format date to Y.m.d H:M:S and convert to epoch time
        date = convert_to_epoch(format_date(i[0]))
        # get instrument
        instrument = map_instrument_type(i[22])
        # keep only novaseq runs if option is set
        to_keep = True
        if keep_novaseq and instrument != 'NovaSeq':
            to_keep = False
        
        if date >= start_date and to_keep:
            if project not in D:
                D[project] = {}
            if run_id not in D[project]:
                D[project][run_id] = {}
            # collect file paths and swids 
            if filename in D[project][run_id]:
                # select files with the most recent workflow id if identical files
                if get_workflow_id(file_path) >= get_workflow_id(D[project][run_id][filename]['filepath']):
                    D[project][run_id][filename]['filepath'] = file_path
                    D[project][run_id][filename]['swid'] = swid
            else:
                D[project][run_id][filename] = {'filepath': file_path, 'swid': swid}
    return D


def get_QC_status_from_nabu(api, file_swid):
    '''
    (str, str) -> (str | None, str)
    
    Returns a tuple with the file QC status and release ticket if file is released
        
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swid (str): File unique identifier
    '''
    
    try:
        response = requests.get(api + '/fileqc/{0}'.format(file_swid), {'accept': 'application/json',})
    except:
        qcstatus, ticket = None, 'NA'
            
    # check response code
    if response.status_code == 200:
        d = response.json()
        if d['fileqcs']:
            assert len(d['fileqcs']) == 1
            qcstatus = d['fileqcs'][0]['qcstatus']
            if 'comment' in d['fileqcs'][0]:
                ticket = d['fileqcs'][0]['comment']
            else:
                ticket = 'NA'
        else:
            qcstatus, ticket = None, 'NA'
    else:
        qcstatus, ticket = None, 'NA'

    return qcstatus, ticket
        


def add_QC_status(api, fastqs):
    '''
    (str, dict) -> dict
    
    Returns a dictionary with file path, file swid, QC status and release ticket organied by project and run
    
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swid (str): File unique identifier
    '''

    for project in fastqs:
        for run in fastqs[project]:
            for filename in fastqs[project][run]:
                # get file QC status
                qcstatus, ticket = get_QC_status_from_nabu(api, fastqs[project][run][filename]['swid'])
                # add qc status and ticket
                fastqs[project][run][filename]['ticket'] = ticket
                fastqs[project][run][filename]['qcstatus'] = qcstatus

    return fastqs


def clean_up_tickets(L):
    '''
    (list) -> list
    
    Returns a list of tickets for a single run or NA if files in the run were
    not released or if QC info could not be retrieved for any of the files in the run
        
    Parameters
    ----------
    - L (list): List of tickets for a single run and/or NA if no release of QC info cannot be extrcted from nabu
    '''
    
    # apply consistent ticket naming scheme by removing URL if present 
    T = list(map(lambda x: os.path.basename(x).upper(), L))
    # remove 'NA' if at least 1 ticket is found.
    # this is because 'NA' is added if the QC status cannot be retrieved from Nabu    
    if any(map(lambda x: x.startswith('GDR'), T)):
        while 'NA' in T:
            T.remove('NA')
    return T


def add_links_to_tickets(L):
    '''
    (list) -> list

    Returns a list of valid Jira URLs or NA if release ticket not defined     

    Parameters
    ----------    
    - L (list): List of tickets for a single run and/or NA if no release of QC info cannot be extrcted from nabu
    '''
    
    T = ['https://jira.oicr.on.ca/browse/{0}'.format(i) if i != 'NA' else i for i in L]
    return T


def write_table(table_file, fastqs):
    '''
    (str, dict) -> None
    
    Parameters
    ----------
    - table_file (str): Path to table file with data release status
    - fastqs (dict): Dictionary with file path, file swid, QC status and release ticket organied by project and run
    '''
    
    newfile = open(table_file, 'w')
    header = ['Project', 'Run', 'File_count', 'Number_files_released', 'Number_files_missing_status', 'Percent_missing_status', 'Released', 'Tickets'] 
    newfile.write('\t'.join(header) + '\n')
     

    # make a sorted list of projects
    projects = sorted(list(fastqs.keys()))
    for project in projects:
        # make a sorted list of runs
        runs = sorted(list(fastqs[project].keys()))
        for run in runs:
            # count number of files
            file_count = len(fastqs[project][run])
            # get the tiket and number of files with PASS status (ie released) + clean up
            # clean up ticket lists: consistent format and remove NA due to missing QC info
            tickets = clean_up_tickets(list(set([fastqs[project][run][filename]['ticket'] for filename in fastqs[project][run]])))
            # add links to tickets
            tickets = add_links_to_tickets(tickets)
            # get release status from ticket list
            if any(map(lambda x: os.path.basename(x).startswith('GDR'), tickets)):
                release_status = 'YES'
            else:
                release_status = 'NO'
            num_released_files =  [fastqs[project][run][filename]['qcstatus'] for filename in fastqs[project][run]].count('PASS')   
            num_missing_status =  len([fastqs[project][run][filename]['qcstatus'] for filename in fastqs[project][run] if fastqs[project][run][filename]['qcstatus'] is None])  
            percent_missing = round(num_missing_status / file_count * 100, 2)
            line = list(map(lambda x: str(x), [project, run, file_count, num_released_files, num_missing_status, percent_missing, release_status, ' '.join(tickets)]))
            newfile.write('\t'.join(line) + '\n')

    newfile.close()


def track_files(args):
    '''
    (list) -> None
    
    Write a tab deliminated table with release information for each run and project for fastqs sequenced at OICR
    by parsing the File Provenance Report and cross-referencing with Nabu
        
    Parameters
    ----------
    
    - fpr (str): Path to File Provenance Report
    - m (int): Number of months prior the current date from which records are considered
    - t (str): Path to table file with data release status 
    - a (str): URL of the Nabu API
    - novaseq (bool): Keep only novaseq runs if activated 
    '''
    
    # dereference FPR
    provenance = os.path.realpath(args.provenance)
    # parse file provenance report
    fastqs = extract_fastqs(provenance, args.interval, args.novaseq)
    # add file QC status from nabu
    fastqs = add_QC_status(args.api, fastqs)
    # write table file
    write_table(args.table, fastqs)


if __name__ == '__main__':
    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'release_tracker.py', description='A tool to track released fastqs')
    parser.add_argument('-fpr', '--provenance', dest='provenance', default='/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    parser.add_argument('-m', '--months', dest='interval', default=12, type=int, help='Number of months prior the current date from which records are considered')
    parser.add_argument('-t', '--table', dest='table', help='Path to table file with data release status')
    parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    parser.add_argument('--novaseq', dest='novaseq', action='store_true', help='Keep only novaseq runs if activated')
    parser.set_defaults(func=track_files)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)