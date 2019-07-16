#!/bin/bash

# DATE: 05-07-19
# TASK: Automated Git tool
# AUTHOR: Matt mcElheron

echo "Jeeves will now get your gits in order"
# First a file is made or changed
# Files are then added to the staging area by using the git add command and passing necessary options
git add .

# Files are then commited to the local repository using the commit command
# -m allows a message to be added
git commit -m "Thanks Jeeves"

# Changes are then sent to GitHub using push
git push
