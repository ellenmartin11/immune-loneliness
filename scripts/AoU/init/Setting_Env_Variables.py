#!/usr/bin/env python
# coding: utf-8

# # Setting and Loading VWB Environment Variables
# 
# This notebook ensures that key workspace-level variables — like the Google Cloud project, workspace buckets, and CDR — are available to every notebook in this Jupyter environment.
# 
# It works by:
# 1. Writing the variables to your `~/.bashrc` file (so they persist across sessions).
# 2. Using a helper function to load them into Python’s `os.environ`.
# 3. Allowing other notebooks to simply import these variables instead of redefining them.
# 
# ---
# 
# ### Variables that will be set
# - `GOOGLE_CLOUD_PROJECT`
# - `WORKSPACE_BUCKET`
# - `WORKSPACE_TEMP_BUCKET`
# - `WORKSPACE_CDR`
# 

# We will use the WB CLI commands to dynamically obtain our environment variables and set them. This is the preferred method as updates to the datasets with each release can occur and overwrite pre-existing data.

# ## Step 1: Read and write environment variables to `.bashrc`
# 
# This step:
# - extract main env variables from wb resource
# - Checks whether `~/.bashrc` exists (creates it if not).
# - Appends export statements for the four variables.
# - If a variable already exists, it replaces its old value with the new one.
# 
# You’ll only need to run this cell once unless you want to change the variable values.

# In[ ]:


import json
import subprocess
import os
import re

# --- Extract variables directly into os.environ ---
workspace = json.loads(subprocess.run(
    ["wb", "workspace", "describe", "--format=json"],
    capture_output=True, text=True, check=True
).stdout)

os.environ["GOOGLE_CLOUD_PROJECT"] = workspace["googleProjectId"]

resources = json.loads(subprocess.run(
    ["wb", "resource", "list", "--format=json"],
    capture_output=True, text=True, check=True
).stdout)

# List to collect all matching CDRs
cdr_matches = []
# Updated Pattern: Starts with either 'C' or 'R'
# Groups: (Prefix), (Year), (Quarter), (Release)
cdr_pattern = re.compile(r"^([CR])(\d{4})Q(\d+)R(\d+)$")

# --- Step 3: Extract workspace resources ---
for r in resources:

 # 1. MASTER BUCKET LOGIC (Outer block)
    if r["resourceType"] == "GCS_BUCKET":
        print(f"Found bucket: id={r['id']}, bucketName={r['bucketName']}")
        
        # Inner conditional logic strictly for filtering different bucket IDs
        if "temporary-workspace-bucket" in r["id"]:
            os.environ["WORKSPACE_TEMP_BUCKET"] = f"gs://{r['bucketName']}"
        elif "workspace-bucket" in r["id"]:
            os.environ["WORKSPACE_BUCKET"] = f"gs://{r['bucketName']}"
        elif r["id"] == "aou-tutorial-notebooks":
            os.environ["bucket_aou_tutorial"] = f"gs://{r['bucketName']}"
            os.environ["bucket_id_aou_tutorial"] = r["id"]
            print("-> Assigned aou-tutorial variables")
        elif r["id"].startswith("rw-migration"):
            os.environ["bucket_migrated"] = f"gs://{r['bucketName']}"
            os.environ["bucket_id_migrated"] = r["id"]
            print(f"-> Assigned migration variables (ID: {r['id']})")
            
    # 2. BQ DATASET LOGIC (Handles C... and R...)
    elif r["resourceType"] in ["BQ_DATASET", "BIGQUERY_DATASET"]:
        dataset_id = r["datasetId"]
        match = cdr_pattern.match(dataset_id)
        
        if match:
            prefix, year, quarter, release = match.groups()
            # The sorting key now includes the prefix to keep C and R distinct if needed,
            # but still prioritizes the latest Year/Quarter/Release.
            version_key = (prefix, int(year), int(quarter), int(release))
            full_path = f"{r['projectId']}.{dataset_id}"
            cdr_matches.append((version_key, full_path))

# --- Finalize the WORKSPACE_CDR ---
if cdr_matches:
    # This sorts by Prefix (C before R), then Year, then Quarter, then Release
    cdr_matches.sort(key=lambda x: x[0])
    latest_r = cdr_matches[-1][1]
    
    os.environ["WORKSPACE_CDR"] = latest_r
    print(f"✅ Successfully identified latest dataset: {os.environ['WORKSPACE_CDR']}")
else:
    os.environ["WORKSPACE_CDR"] = ""
    print("⚠️ No C... or R... datasets found.")

# --- Save to .bashrc ---
# --- Expanded tracking array containing all variables ---
vars_to_save = [
    "GOOGLE_CLOUD_PROJECT", 
    "WORKSPACE_BUCKET", 
    "WORKSPACE_TEMP_BUCKET", 
    "WORKSPACE_CDR",
    "bucket_aou_tutorial",     # Path (gs://...)
    "bucket_id_aou_tutorial",  # ID only (aou-tutorial-notebooks)
    "bucket_migrated",         # Path (gs://...)
    "bucket_id_migrated"       # ID only (rw-migration-xxxx)
]

# --- Print validation check ---
print("\nVariables extracted:")
for var in vars_to_save:
    value = os.environ.get(var)
    print(f"{var}: {value if value else 'NOT FOUND'}")


# --- Save to .bashrc ---
bashrc_path = os.path.expanduser("~/.bashrc")

if not os.path.exists(bashrc_path):
    print(f"Creating {bashrc_path}...")
    with open(bashrc_path, "w") as f:
        f.write("# Created by Verily setup script\n")

# Append everything safely to the end of the file
with open(bashrc_path, "a") as f:
    f.write("\n# Verily Workbench variables\n")
    for var in vars_to_save:
        value = os.environ.get(var)
        if value:
            f.write(f'export {var}="{value}"\n')

print(f"\n✅ Saved to {bashrc_path}")


# In[ ]:


for r in resources:
    print(r["id"], r["resourceType"])
    if r["resourceType"] == "GCS_BUCKET" and "workspace-bucket" in r["id"]:
        os.environ["WORKSPACE_BUCKET"] = f"gs://{r['bucketName']}"
    elif r["resourceType"] == "GCS_BUCKET" and "temporary" in r["id"]:
        os.environ["WORKSPACE_TEMP_BUCKET"] = f"gs://{r['bucketName']}"


# ## Step 2: Load Environment Variables in Python and test
# 
# The following function reads your `~/.bashrc` file and sets each variable in Python’s runtime environment (`os.environ`), so you can use them directly in your code.
# 
# Any other notebook can reuse this function to reload environment variables at the start.

# In[ ]:


import os

# In your notebook, just run this:
with open(os.path.expanduser("~/.bashrc"), 'r') as f:
    for line in f:
        if line.strip().startswith('export '):
            parts = line.strip().replace('export ', '').split('=', 1)
            if len(parts) == 2:
                var_name = parts[0].strip()
                var_value = parts[1].strip().strip("'\"")
                
                # SKIP PATH completely!
                if var_name == 'PATH':
                    continue  # Skip this line
                    
                os.environ[var_name] = var_value

# Now use them
print(f"GOOGLE_CLOUD_PROJECT = {os.environ.get('GOOGLE_CLOUD_PROJECT')}")
print(f"WORKSPACE_BUCKET = {os.environ.get('WORKSPACE_BUCKET')}")


# **Confirm that the environment variables are correctly set.**

# In[ ]:


#!echo $WORKSPACE_CDR


# In[ ]:


#!echo $WORKSPACE_TEMP_BUCKET


# In[ ]:


# Now you can access them directly:
print(f"WORKSPACE_CDR = {os.environ.get('WORKSPACE_CDR')}")
print(f"WORKSPACE_BUCKET = {os.environ.get('WORKSPACE_BUCKET')}")
print(f"GOOGLE_PROJECT = {os.environ.get('GOOGLE_PROJECT')}")
