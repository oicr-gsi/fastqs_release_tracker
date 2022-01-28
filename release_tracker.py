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
    
    Returns the workflow id from the file path
           
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
    
    
    '''
    
    try:
        records = subprocess.check_output('zcat {0} | grep {1}'.format(provenance, workflow), shell=True).decode('utf-8').rstrip().split('\n')
    except:
        records = []
    finally:
        return records


def extract_fastqs(provenance, time_interval):
    '''
    (str, int) -> (dict)
  
    Returns xxxxx
            
    Parameters
    ----------
    - provenance (str): Path to File Provenance Report
    - time_interval (int): Number of months prior the current date from which records are considered
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
            
        if date >= start_date:
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



# def get_project_QC_status(api, project, workflow):
#     '''
    
    

#     Parameters
#     ----------
#     api : TYPE
#         DESCRIPTION.
#     project : TYPE
#         DESCRIPTION.
#     workflow : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     '''


#     headers = {
#     'accept': 'application/json',
# }

# params = (
#     ('project', 'HCCCFD'),
#     ('workflow', 'blc2fastq'),
#     ('qcstatus', 'PASS'),
# )

# response = requests.get('http://gsi-dcc.oicr.on.ca:3000/fileqcs', headers=headers, params=params)




#     try:
#         response = requests.get(api + '/fileqcs?project={0}&workflow={1}'.format(project, workflow))
#     except:
#         raise ValueError('Cannot retrieve QC status from Nabu')
    
#     # check response code
#     if response.status_code == 200:
#         d = response.json()
#         assert len(d['fileqcs']) == 1
#         filepath = d['fileqcs'][0]['filepath']
#         fileswid = d['fileqcs'][0]['fileswid']
#         project = d['fileqcs'][0]['project']
#         qcstatus = d['fileqcs'][0]['qcstatus']
#         if 'comment' in d['fileqcs'][0]:
#             ticket = d['fileqcs'][0]['comment']
#         else:
#             ticket = 'NA'
#     else:
#         raise ValueError('Could not get file QC status from nabu for {0}. Nabu response code: {1}'.format(file_swid, response.status_code))

#     return qcstatus, ticket, filepath, fileswid, project
    






def get_QC_status_from_nabu(api, file_swid):
    '''
    (str, str, str, str, str | None) -> None
    
    
    XXXXXXXXXXXXXXXXXXXXXXx
    
    
    
    
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
        assert len(d['fileqcs']) == 1
        qcstatus = d['fileqcs'][0]['qcstatus']
        if 'comment' in d['fileqcs'][0]:
            ticket = d['fileqcs'][0]['comment']
        else:
            ticket = 'NA'
    else:
        qcstatus, ticket = None, 'NA'

    return qcstatus, ticket
        


def add_QC_status(api, fastqs):
    '''
    
    
    

    Returns
    -------
    None.

    '''

    for project in fastqs:
        print(project)
        for run in fastqs[project]:
            for filename in fastqs[project][run]:
                # get file QC status
                qcstatus, ticket = get_QC_status_from_nabu(api, fastqs[project][run][filename]['swid'])
                # add qc status and ticket
                fastqs[project][run][filename]['ticket'] = ticket
                fastqs[project][run][filename]['qcstatus'] = qcstatus

    return fastqs



def write_table(table_file, fastqs):
    '''
    
    
    '''
    
    
    newfile = open(table_file, 'w')
    header = ['project', 'run', 'file_count', 'number_files_released', 'number_files_missing_status', 'percent_missing_status', 'tickets'] 
    newfile.write('\t'.join(header) + '\n')
    # project run files_num ticket num_files_reeleased 
    

    # make a sorted list of projects
    projects = sorted(list(fastqs.keys()))

    for project in projects:
        # make a sorted list of runs
        runs = sorted(list(fastqs[project].keys()))
        for run in runs:
            # count number of files
            file_count = len(fastqs[project][run])
            # get the tiket and number of files with PASS status (ie released)
            tickets = list(set([fastqs[project][run][filename]['ticket'] for filename in fastqs[project][run]]))
            num_released_files =  [fastqs[project][run][filename]['qcstatus'] for filename in fastqs[project][run]].count('PASS')   
            num_missing_status =  len([fastqs[project][run][filename]['qcstatus'] for filename in fastqs[project][run] if fastqs[project][run][filename]['qcstatus'] is None])  
            percent_missing = round(num_missing_status / file_count * 100, 4)
            line = list(map(lambda x: str(x), [project, run, file_count, num_released_files, num_missing_status, percent_missing, ';'.join(tickets)]))
            newfile.write('\t'.join(line) + '\n')

    newfile.close()


def track_files(args):
    
    # dereference FPR
    provenance = os.path.realpath(args.provenance)
    fastqs = extract_fastqs(provenance, args.interval)
    fastqs = add_QC_status(args.api, fastqs)
    write_table(args.table, fastqs)




    
    
if __name__ == '__main__':

    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'release_tracker.py', description='A tool to track released fastqs')
    parser.add_argument('-fpr', '--provenance', dest='provenance', default='/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    parser.add_argument('-m', '--months', dest='interval', default=12, type=int, help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    
    
    
    parser.add_argument('-t', '--table', dest='table', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    
    
    parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
                        
                        
                        
    
    
    
    parser.set_defaults(func=track_files)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)