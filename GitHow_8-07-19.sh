#!/bin/bash

# DATE: 05-07-19
# TASK: Usage of GitHub
# AUTHOR: Matt mcElheron

# The user and email was configurated to the system
git config --global user.email "m.m@ucd"
git config --global user.name "m.m@ucd"

# clone is used to add a dir to the used system
git clone link_from_GitHub

# First a file is made or changed
# Files are then added to the staging area by using the git add command and passing necessary options
git add .

# Files are then commited to the local repository using the commit command
# -m allows a message to be added
git commit -m "hello future me"





