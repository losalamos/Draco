import requests
import json
import re
import os

#Bypass LANL proxy
os.environ['no_proxy'] = '*'

#This is the gitlab private access token (different from GitHub token)
private_token = ''

#GitLab project id as a string, Can retrieve under project name on project page
project_id = '302'

#Gitlab merged url. Need to authenticate in URL due to the limited function of the Gitlab web API.
api_url = 'https://gitlab.lanl.gov/api/v4/projects/{0}/merge_requests?state=merged&private_token={1}'.format(project_id, private_token)

result = requests.get(api_url)

print(result)

result_json = result.json()

if(private_token == '' or project_id == ''):
    print("Either your private token or project id is empty.")
    exit()

#Collects and prints information on the first JSON page
for entry in result_json:
    if(entry['merged_at']):
        #Prevents empty string error
        if(entry["merged_at"] != None):
            year, month, day = entry["merged_at"][:10].split('-')
            year, month, day = int(year), int(month), int(day)
        if (year >= 2019 and ((month >= 2 and day >= 26) or month >= 3)):
            print("{0} {1} {2}".format(entry["iid"], entry["title"], entry["web_url"]))

#Collects and prints information on the rest of the JSON pages
while 'next' in result.links.keys():
    result=requests.get(result.links['next']['url'])
    next_result_json = result.json()
    if(entry['merged_at']):
        for entry in next_result_json:
            #Prevents empty string error
            if(entry["merged_at"] != None):
                year, month, day = entry["merged_at"][:10].split('-')
                year, month, day = int(year), int(month), int(day)
            if (year >= 2019 and ((month >= 2 and day >= 26) or month >= 3)):
                print("{0} {1} {2}".format(entry["iid"], entry["title"], entry["web_url"]))
