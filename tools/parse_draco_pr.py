import requests
import json
import re
import os

#Bypass LANL proxy
#os.environ['no_proxy'] = '*'

#This is your github login username
github_username = ''

#This is your github personal access token (Different from GitLab)
github_token = ''

#GitHub closed url.
api_url = 'https://api.github.com/repos/lanl/Draco/pulls?state=closed'

result = requests.get(api_url, auth=(github_username, github_token))

print(result)

result_json = result.json()

if(github_username == '' or github_token == ''):
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
            print("{0} {1} {2}".format(entry["number"], entry["title"], entry["url"]))

#Collects and prints information on the rest of the JSON pages
while 'next' in result.links.keys():
    result=requests.get(result.links['next']['url'], auth=(github_username, github_token))
    next_result_json = result.json()
    if(entry['merged_at']):
        for entry in next_result_json:
            #Prevents empty string error
            if(entry["merged_at"] != None):
                year, month, day = entry["merged_at"][:10].split('-')
                year, month, day = int(year), int(month), int(day)
            if (year >= 2019 and ((month >= 2 and day >= 26) or month >= 3)):
                print("{0} {1} {2}".format(entry["number"], entry["title"], entry["url"]))
