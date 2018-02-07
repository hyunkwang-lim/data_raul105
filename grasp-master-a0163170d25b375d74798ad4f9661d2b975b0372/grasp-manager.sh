#!/bin/bash

# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

### DEVELOPER GUIDE: ###

## How to add a parameter ##

# 1. Add a default value in the block "Definition of default values of parameters"
# 2. Add a documentation in grasp-manayer.yml.dist
# 3. Modify prepare_environment function to read new parameter. It is just to add a block for get that value from grasp-manager.yml
# 4. Use the new parameter as you wish. The developer can access to the value of that parameter in all actions of the script to customize a behavior. 

## How to add a module ##

# 1. Add a default value in the block "Definition of default values of parameters"
# 2. Add a documentation in grasp-manayer.yml.dist
# 3. Add a variable with default_git_server_module_{name}
# 4. Add a default value if the parameter is not present in setting file. (Recommended false) 
# 5. Modify prepare_environment function to read new parameter. It is just to add a block for get that value from grasp-manager.yml
# 6. Add them to repos_to_install_* arrays in order to use it in all actions automatically

export LANGUAGE='en_US.UTF-8 git'

# grasp-manager echo function. It will print in a pretty way
function gmecho_title () {
    echo -e $'\e[38;5;208m'"#### $1 ####"$'\e[0m'
} 

function gmecho () {
    echo -e $'\e[38;5;208m'"$1"$'\e[0m'"$2"
}

if [ ! -f "grasp-manager.yml.dist" ]
then
    gmecho_title "ERROR: grasp-manager script has to be run in a valid grasp root folder."
    gmecho_title "Please, move your terminal to correct place and try again."
    exit -1
fi

# Preparing environment
core_drivers="sdata"
core_transformers="none segment_imagedat" 
core_segment_functions="ascii classic classic_plot none"
core_current_functions="none"
core_tile_functions="ascii none"
core_kernels="KERNELS_BASE"
core_constants_sets="generic valgrind"

compilation_environment=""

default_git_server=git@code.grasp-open.com
default_git_driver_prefix="open/driver-"
default_git_transformer_prefix="open/transformer-"
default_git_segment_function_prefix="open/segment_function-"
default_git_current_function_prefix="open/current_function-"
default_git_tile_function_prefix="open/tile_function-"
default_git_kernels_prefix="open/kernels-"
default_git_constants_set_prefix="open/constants-set-"

default_git_server_module_dependencies=git@code.grasp-open.com:grasp-internal/module-dependencies.git
default_git_server_module_gpgpu=git@code.grasp-open.com:grasp-internal/module-gpgpu.git
default_git_server_module_kernels_modifier=git@code.grasp-open.com:grasp-internal/module-kernels_modifier.git
default_git_server_module_internal_examples=git@code.grasp-open.com:grasp-internal/module-internal_examples.git
default_git_server_module_models=git@code.grasp-open.com:grasp-sas/module-models.git

timestamp=$(date "+%Y%m%d_%H%M%S")

base_path=$(pwd)

home_manager=$base_path"/home/grasp-manager/"

mkdir -p $home_manager 2> /dev/null

bk_folder_path=$home_manager/backup/

last_backup=$(ls -1 $bk_folder_path 2> /dev/null | sort -r | head -n 1)

if [ "$last_backup" != "" ]
then
    last_environment=$(cat $bk_folder_path$last_backup/log.txt |head -n1 |cut -d' ' -f2)
else
    last_environment="no_defined"
fi

bk_path=$bk_folder_path$timestamp/

logfilename=log.txt

logfile=$bk_path/$logfilename

repos_path=$home_manager/repos/

password=""

sudo=""

log_git_token="#restore_info"

# Definition of default values of parameters
compilation_compile_dependencies="true"
compilation_use_sudo="true"
compilation_constants_set="generic"
compilation_build="Release"
compilation_mpi="off"
compilation_f90="gfortran"
compilation_prefix="/usr/local"
compilation_build_dir="build"
compilation_cc="cc"
compilation_cxx="c++"
compilation_cmake_flags=""
modules_dependencies="false"
modules_gpgpu="false"
modules_kernels_modifier="false"
modules_internal_examples="false"
modules_models="false"
retrieval_modules_enabled=""

compilation_flags="" # It is set in prepare_environment function

last_version_tag=$(git tag | grep ^v | sort -t. -k1.2,1nr -k2,2nr -k3,3nr | head -n1) #last version number

# removing old backups. We keep maximum 40
total_backups=$(ls -1 $bk_folder_path 2> /dev/null | wc -l)
total_backups=$(echo $total_backups)
let nbackups_to_remove=$total_backups-100
backups_to_remove=$(ls -1 $bk_folder_path 2> /dev/null | sort -r| tail -n $nbackups_to_remove)
if [ "$nbackups_to_remove" -gt "0" ]
then
    for backup_to_remove in $backups_to_remove
    do
        rm -rf $bk_folder_path$backup_to_remove
    done
fi

# Cleaning homa_manager. We remove everything because with following command
# folders can not be removed. Rest of files can be temporal files
cd $home_manager
rm * 2> /dev/null
cd - > /dev/null

# get current version from git. It can be a branch name or
# a sentences type * (HEAD detached at c93f004) that we can 
current_version=$(git branch | grep '*' |cut -d' ' -f2)
if [ "${current_version:0:1}" = "(" ] 
then
    current_version=$(git branch | grep '*')
    current_version=$(echo ${current_version:19}|cut -d')' -f1)
fi

if [ ! -e grasp-manager.yml ]
then
    cp grasp-manager.yml.dist grasp-manager.yml
fi

# Definition of general functions

# Print help information
function help {
    gmecho_title "GRASP MANAGER"
    echo "This script helps you to manage grasp code: updating, compiling, changing version..."
    echo ""
    echo "Usage:"
    echo "  ./grasp-manager [action] [action-arguments] [environment]"
    echo ""
    echo "When a enviroment is not defined last environment is used for configuration definition"
    echo "and current installed extensions as required extensions" 
    echo ""
    echo ""
    echo "List of available actions:"
    echo "--------------------------"
    echo ""
    echo "For users:"
    gmecho "  update-grasp [environment]" ": This action backed up your changes in code and download, compile and install latest stable version of the code"
    gmecho "  update-grasp-to [version] [environment]" ": This action backed up your changes in code and download, compile and install version required as first argument"
    gmecho "  update-grasp-to-dev [environment]" ": This action backed up your changes in code and download, compile and install lastest version of the code (UNSTABLE)"    
    gmecho "  purge [all-backups|unused-repos|all]" ": This action will remove grasp-manager cache files. User can choose if he want to remove all-backups, only cached repos not present in grasp-manager.yml file or both options"
    gmecho "  install-extension [extension-type] [extension-name] [extension-url]" ": This action install an single external repository. extension-type can be 'driver', 'transformer', 'segment_function', 'current_function', 'tile_function', 'kernels', 'constants_set' or 'module'. extension-url can be omitted and default value will be used"
    gmecho "  uninstall-extension [extension-type] [extension-name]" ": This action install an single external repository. extension-type can be 'driver', 'transformer', 'segment_function', 'current_function', 'tile_function', 'kernels', 'constants_set' or 'module'"
    echo ""
    echo "For developers:"
    gmecho "  make [environment]" ": Just compile the code using 'environment' constants. If environment is not defined last used will be used (or default if it does not exist)"    
    gmecho "  pull [environment]" ": It synchronizes the code of 'environment'. If environment is not defined installed extensions will be updated"
    gmecho "  checkout [commit, tag or branch name] [environment]" ": This action perform a checkout to a specfic commit moving all extensions there. If environment is not defined installed extensions will be updated"
    gmecho "  clean" ": Just clean compiled code. It is an alias to 'make clean' command"        
    gmecho "  apply [backup-path]" ": This action comes back to a specific version backed up and change environment for the used in that version."
    gmecho "  rollback" ": Undo last action applying last backup. As 'apply' action but over last backup"
    gmecho "  commit [message]" ": It performs same commit (adding all changes and pushing them) in core code and in all installed extensions"
    gmecho "  status" ": It shows git status command for all repositories (core and installed extensions)"
    gmecho "  tag [tag-name]" ": It applies a tag to all installed extension and core in its current GIT commit"
    gmecho "  merge [branch]" ": It merges 'branch' in core code and in all installed extensions"
    gmecho "  branch [branch]" ": It creates a branch called 'branch' in core code and in all installed extensions"
    gmecho "  push" ": It perform a git push in core code and all installed extensions"
    gmecho "  create-extension [extension-type] [extension-name] [extension-url]" ": It creates a new extension. extension-type can be 'driver', 'transformer', 'segment_function', 'current_function' or 'tile_function'. If extension-url is not provided a default value will be used"
    echo ""
}

function backup {    
    mkdir -p $bk_path 2> /dev/null

    # Saving current changes and undo them
    changes=$(git status -s)
    if [ "$changes" != "" ]
    then
        git add -A .
        git stash
        git stash show stash@{0} -p --binary --color=never > $bk_path/core.patch
        git stash drop > /dev/null
        gmecho_title "changes have been backed up in $bk_path/core.patch"
        gmecho_title "You can restore them with 'git stash apply' or using the patch file"
    fi

    for i in ${!installed_repo_names[@]}
    do
        cd ${installed_repo_paths[$i]}
        remote_url=$(git config --get remote.origin.url)
        guess_path=${remote_url/:/\/}
        cd - > /dev/null

        if [ -d "${installed_repo_paths[$i]}/.git" -a -d "$repos_path/$guess_path" ]
        then
            cd ${installed_repo_paths[$i]}
            changes=$(git status -s)
            if [ "$changes" != "" ]
            then
                git add -A .
                git stash
                git stash show stash@{0} -p --binary --color=never > $bk_path/${installed_repo_types[$i]}-${installed_repo_names[$i]}.patch
                git stash drop > /dev/null
                echo "changes in ${installed_repo_names[$i]} have been backed up in $bk_path/${installed_repo_types[$i]}-${installed_repo_names[$i]}.patch"
                echo "You can restore them with 'git stash apply' or using the patch file"
            fi
            cd - > /dev/null
        fi
    done
}

function make_dependencies {
    if [ "$compilation_compile_dependencies" = "true" ]
    then
        if [ "$compilation_use_sudo" = "true" ]
        then
            apply_sudo=$sudo
        else
            apply_sudo=""
        fi
        bash -c "$apply_sudo make -f Makefile dep $compilation_flags"
    fi
}

function make_clean {
    if [ "$password" != "" ]
    then
        apply_sudo=$sudo
    else
        apply_sudo=""
    fi
    bash -c "$apply_sudo make -f Makefile clean"
}

function make_compile {
    make_clean
    bash -c "make -f Makefile $compilation_flags"
}

function make_install {
    if [ "$compilation_use_sudo" = "true" ]
    then
        apply_sudo=$sudo
    else
        apply_sudo=""
    fi
    bash -c "$apply_sudo make -f Makefile install $compilation_flags"
}

function manager_purge () {
    if [ "$1" != "all-backups" -a "$1" != "all-repos" -a "$1" != "unused-repos" -a "$1" != "all" -a "$1" != "manual-installed" -a "$1" != "all-installed-repos" ]
    then
        echo "ERROR: in manager_purge function"
        exit 0
    fi

    # Cleaning all backups
    if [ "$1" = "all-backups" -o "$1" = "all" ]
    then
        cd $bk_folder_path
        rm -rf *
        cd - > /dev/null
    fi

    # Secret function. Not documented. Cleaning all repos
    if [ "$1" = "all-repos" ]
    then
        cd $repos_path
        rm -rf *
        cd - > /dev/null
    fi

    # Cleaning unused repos
    if [ "$1" = "unused-repos" -o "$1" = "all" ]
    then
        cd $repos_path
        repos_list=$(find ./ -maxdepth 3 -mindepth 1 -type d | cut -c 4-999999 | sed "s/\//:/")
        for repo in $repos_list
        do
            grep -w $repo ../../../grasp-manager.yml > /dev/null
            if [ "$?" != "0" ] 
            then
                rm -rf $repo
            fi
        done
        cd - > /dev/null
    fi

    # It is not documented because it should not be used
    if [ "$1" = "manual-installed" ]
    then
        prepare_environment
        for i in ${!installed_repo_names[@]}
        do
            cd ${installed_repo_paths[$i]}
            remote_url=$(git config --get remote.origin.url)
            guess_path=${remote_url/:/\/}
            cd - > /dev/null

            if [ -d "${installed_repo_paths[$i]}/.git" -a -d "$repos_path/$guess_path" ]
            then
                echo -ne "" > /dev/null # Do nothing
            else
                gmecho_title "Removing ${installed_repo_types[$i]} ${installed_repo_names[$i]} (Moving to backups)"
                mv ${installed_repo_paths[$i]} $bk_path/${installed_repo_types[$i]}-${installed_repo_names[$i]}
            fi
        done
    fi

    # It is not documented because it should not be used
    if [ "$1" = "all-installed-repos" ]
    then
        prepare_environment
        for i in ${!installed_repo_names[@]}
        do
            rm -rf ${installed_repo_paths[$i]}
        done
    fi
}

# This function will install an external extension
# $1 type of extension
# $2 name of the extension
# $3 If git url has to be provided  Example manager_install_extension kernel parasol git@code.grasp-open... 
function manager_install_extension () {
    if [ "$1" != "driver" -a "$1" != "transformer" -a "$1" != "segment_function" -a "$1" != "current_function" -a "$1" != "tile_function" -a "$1" != "module" -a "$1" != "kernels" -a "$1" != "constants_set" ]
    then
        echo "ERROR: Extension type has to be 'driver' or 'transformer' or 'segment_function' or 'current_function or 'tile_function' or 'module' or 'kernels' or 'constants_set'"
        exit 0
    fi

    if [ "$2" != "${2/-/}" ]
    then
        echo "ERROR: Argument 2 (extension name) can not contain dash symbol (-). Subtitute it by an underscore (_)"
        exit -1
    fi

    # Extracting the information from the command line
    extension_type=$1
    extension_name=${2/:/}
    extension_url=$3
    if [ "$extension_url" = "" ]
    then
        extension_url="~"
    fi

    # Preparing the temporal file
    tmp_yaml_file=$home_manager"/tmp-grasp-manager.yml"

    case "$extension_type" in
        "driver")
            echo -e "$last_environment:\n    extensions:\n        input:\n            drivers:\n                $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "transformer")
            echo -e "$last_environment:\n    extensions:\n        input:\n            transformers:\n                $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "segment_function")
            echo -e "$last_environment:\n    extensions:\n        output:\n            segment_functions:\n                $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "current_function")
            echo -e "$last_environment:\n    extensions:\n        output:\n            current_functions:\n                $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "tile_function")
            echo -e "$last_environment:\n    extensions:\n        output:\n            tile_functions:\n                $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "kernels")
            echo -e "$last_environment:\n    kernels:\n        $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "constants_set")
            echo -e "$last_environment:\n    constants_sets:\n        $extension_name: $extension_url" > $tmp_yaml_file
            ;;
        "module")
            echo -e "$last_environment:\n    modules:\n        $extension_name: true" > $tmp_yaml_file
            ;;
    esac

    # Call to install action
    prepare_environment "$last_environment" "$tmp_yaml_file"

    rm $tmp_yaml_file

    # installing extension
    install_extensions

    git_branch=$(git rev-parse --abbrev-ref HEAD)

    if [ "${repos_to_install_paths[0]}" != "" ]
    then
        if [ "$git_branch" != "HEAD" ]
        then
            # set extension to same branch if it existe
            gmecho_title "Changing version and updating ${repos_to_install_types[0]} ${repos_to_install_names[0]}"
            cd ${repos_to_install_paths[0]}
            extension_tag=$(select_extension_version $git_branch)
            git checkout $extension_tag | grep -v "Branch dev set up to track remote branch dev from origin." | grep -v "Your branch is up-to-date with 'origin/dev'." 
            git pull
            cd - > /dev/null
        fi
    fi
}

# This function will uninstall an external extension
# $1 type of extension
# $2 name of the extension.
function manager_uninstall_extension () {
    if [ "$1" != "driver" -a "$1" != "transformer" -a "$1" != "segment_function" -a "$1" != "current_function" -a "$1" != "tile_function" -a "$1" != "module" -a "$1" != "kernels" -a "$1" != "constants_set" ]
    then
        echo "ERROR: Extension type has to be 'driver' or 'transformer' or 'segment_function' or 'current_function or 'tile_function' or 'module' or 'kernels' or 'constants_set'"
        exit 0
    fi

    if [ "$2" != "${2/-/}" ]
    then
        echo "ERROR: Argument 2 (extension name) can not contain dash symbol (-). Subtitute it by an underscore (_)"
        exit -1
    fi

    # Extracting the information from the command line
    extension_type=$1
    extension_name=${2/:/}

    # Preparing the temporal file
    tmp_yaml_file=$home_manager"/tmp-grasp-manager.yml"

    echo -e "$last_environment: ~" > $tmp_yaml_file
    
    # Call to install action
    prepare_environment "$last_environment" "$tmp_yaml_file"

    rm $tmp_yaml_file

    # We just uninstall required repo
    for i in ${!repos_to_uninstall_names[@]}
    do
        if [ "${repos_to_uninstall_types[$i]}" != "$extension_type" -o "$extension_name" != "${repos_to_uninstall_names[$i]}"  ]
        then
            unset repos_to_uninstall_names[$i]
            unset repos_to_uninstall_types[$i]
            unset repos_to_uninstall_paths[$i]
        fi
    done

    uninstall_extensions
}


function uninstall_extensions {
    if [ "$using_last_environment" = "no" ]
    then
        for i in ${!repos_to_uninstall_names[@]}
        do
            # To know if it is repos folder I get where we should find it
            cd ${repos_to_uninstall_paths[$i]}
            remote_url=$(git config --get remote.origin.url)
            guess_path=${remote_url/:/\/}
            cd - > /dev/null

            # Checking if repo was installed manually
            # A repo is installed manually if it has not .git folder (it is not a repo)
            # or it has git repo but this repo is not present in grasp-manager repos folder
            if [ -d "${repos_to_uninstall_paths[$i]}/.git" -a -d "$repos_path/$guess_path" ] 
            then
                gmecho_title "Uninstalling ${repos_to_uninstall_types[$i]} ${repos_to_uninstall_names[$i]}"
                # It is moved to repos_path but in steps in order to avoid as maximum as possible problems if the user interrupt the script
                # Additionally I prevent CTRL-C keys in that action
                trap 'mv $repos_path/${guess_path}__tmp ${repos_to_uninstall_paths[$i]} && exit -1' 2
                mv ${repos_to_uninstall_paths[$i]} $repos_path/${guess_path}__tmp
                trap 'rm -rf $repos_path/$guess_path && mv $repos_path/${guess_path}__tmp $repos_path/$guess_path && exit -1' 2
                rm -rf $repos_path/$guess_path
                trap 'mv $repos_path/${guess_path}__tmp $repos_path/$guess_path && exit -1' 2
                mv $repos_path/${guess_path}__tmp $repos_path/$guess_path
                trap - 2
            else
                gmecho_title "WARNING: Skip uninstallation of ${repos_to_uninstall_types[$i]} ${repos_to_uninstall_names[$i]} because it was installed manually"
            fi
        done
    fi
}

# It install all extension
function install_extensions {
    if [ "$using_last_environment" = "no" ]
    then
        for i in ${!repos_to_install_names[@]}
        do
            gmecho_title "Installing ${repos_to_install_types[$i]} ${repos_to_install_names[$i]}"

            repo_path=${repos_to_install_urls[$i]/:/\/}

           # If the repo does not exist, we download it
           if [ ! -d "$repos_path/$repo_path" ]
           then
                status_git=1
                tries=0
                while [ "$status_git" != "0" -a "$tries" -lt "10" ]
                do
                    if [ "$tries" -gt "0" ]
                    then
                        gmecho_title "WARNING: Repository could not be cloned. Trying again ($tries tries)"
                    fi 
                    git clone ${repos_to_install_urls[$i]} $repos_path/$repo_path
                    status_git=$?
                    let tries=$tries+1
                done 
                if [ "$tries" -ge "10" ]
                then
                    gmecho_title "ERROR: Repository could not be cloned after $tries tries. grasp-manager will stop"
                    exit -1
                fi
           fi
           trap 'rm -rf ${repos_to_install_paths[$i]} && exit -1' 2
           cp -a $repos_path/$repo_path ${repos_to_install_paths[$i]}
           trap - 2
        done
    fi
}

# Download new code from grasp core and extensions. It updates remote branches
# This function update last_version variable
function git_fetch {
    gmecho_title "Downloading new core code"
    git fetch --all | grep -v "Fetching origin"

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Downloading new code from ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            git fetch --all | grep -v "Fetching origin"
            cd - > /dev/null
        else
            gmecho_title "WARNING: The ${installed_repo_types[$i]} ${installed_repo_names[$i]} is not a git repository and it can not be updated"
        fi
    done

    last_version_tag=$(git tag | grep ^v | sort -t. -k1.2,1nr -k2,2nr -k3,3nr | head -n1) #last version number
}

# Merge a branch in current
function git_merge () {
    gmecho_title "Merging core code"
    git merge "$1"

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Merging code from ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            git merge "$1"
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}

# Merge a branch in current
function git_branch () {
    gmecho_title "Creating branch in core reposiroty"
    git branch "$1"
    git checkout "$1"
    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Creating branch in ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            git branch "$1"
            git checkout "$1"
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}


# Merge remote branch in local branch if it is possible
# Update local branch from remote
# TODO: Currently merge is used to synch with remote. In the future could be possible to extend this feature to rebase too.
function git_update_from_remote {
    gmecho_title "Merging new core code from remote"
    remote=$(git rev-parse --abbrev-ref --symbolic-full-name @{u} 2> /dev/null)
    if [ "$?" = "0" ] # there is remote
    then
        git merge $remote
    fi

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Merging new code from remote of ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            remote=$(git rev-parse --abbrev-ref --symbolic-full-name @{u} 2> /dev/null)
            if [ "$?" = "0" ] # there is remote
            then
                git merge $remote
            fi
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}

# Adding code (git add) in all repositories
function git_add {
    gmecho_title "Adding changes in core repository"
    git add -A .
    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Adding changes in ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            git add -A .
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}

# Committing changes in all repositories
function git_commit () {
    gmecho_title "Committing changes in core repository"
    git commit -m "$1"

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Committing changes in ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            git commit -m "$1"
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}

# Tagging all repositories
function git_tag () {
    gmecho_title "Tagging core repository"
    tagversion=$1
    git tag -a $tagversion -m "Created tag $tagversion"
    git push origin $tagversion

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Tagging ${installed_repo_types[$i]} ${installed_repo_names[$i]} repository"
            cd ${installed_repo_paths[$i]}
            tagversion="$1.0"
            git tag -a $tagversion -m "Created tag $tagversion"
            git push origin $tagversion
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}


# Call git status for all repositories
function git_status () {
    if [ ${#repos_to_install_names[@]} != "0" -o ${#repos_to_uninstall_names[@]} != "0" ]
    then
        gmecho "WARNING: repos are not synchronized. Try to run 'pull' action to synch all repos"
    fi
    if [ ${#repos_to_install_names[@]} != "0" ]
    then
        repos=""
        for i in ${!repos_to_install_names[@]}
        do
            repos=$repos" "${repos_to_install_names[$i]} 
        done        
        gmecho "  There are repositories to be installed ($repos )"
    fi
    if [ ${#repos_to_uninstall_names[@]} != "0" ]
    then
        repos=""
        for i in ${!repos_to_uninstall_names[@]}
        do
            repos=$repos" "${repos_to_uninstall_names[$i]} 
        done        
        gmecho "  There are repositories to be uninstalled ($repos )"
    fi
    if [ ${#repos_to_install_names[@]} != "0" -o ${#repos_to_uninstall_names[@]} != "0" ]
    then
        echo ""
    fi

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Git status for ${installed_repo_types[$i]} function ${installed_repo_names[$i]} repository"
            cd ${installed_repo_paths[$i]}
            git status
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
    
    # For core repository
    gmecho_title "Git status for core repository"
    git status
}

# Pushing changes in all repositories
function git_push () {
    gmecho_title "Pushing changes in core repository"
    git push

    for i in ${!installed_repo_names[@]}
    do
        if [ -d "${installed_repo_paths[$i]}/.git" ]
        then
            gmecho_title "Pushing changes in ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
            cd ${installed_repo_paths[$i]}
            git push
            cd - > /dev/null
        else
            gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
        fi
    done
}

function git_checkout () {
    tag=$1

    gmecho_title "Changing version of core code"
    git checkout $tag 
    if [ $? = 0 ]
    then
        for i in ${!installed_repo_names[@]}
        do
            if [ -d "${installed_repo_paths[$i]}/.git" ]
            then
                gmecho_title "Changing version of ${installed_repo_types[$i]} ${installed_repo_names[$i]}"
                cd ${installed_repo_paths[$i]}
                extension_tag=$(select_extension_version $tag)
                git checkout $extension_tag | grep -v "Branch dev set up to track remote branch dev from origin." | grep -v "Your branch is up-to-date with 'origin/dev'." 
                cd - > /dev/null
            else
                gmecho_title "WARNING: Skipping ${installed_repo_types[$i]} ${installed_repo_names[$i]}. It is not a git repository"
            fi
        done  
    else
        gmecho_title "Error changing version of the code"
    fi
}


# This function write initial status for repository in a log
# this function is called by prepare environment
function write_log {
    mkdir -p $bk_path 2> /dev/null
    touch $logfile

    echo "ENVIRONMENT: ${compilation_environment}" >> $logfile
    git_branch=$(git rev-parse --abbrev-ref HEAD)
    git_commit=$(git rev-parse HEAD)

    echo "#core_restore_info GRASP-core $git_branch $git_commit" >> $logfile


    for i in ${!installed_repo_names[@]}
    do
        # To know if it is repos folder I get where we should find it
        cd ${installed_repo_paths[$i]}
            path_absolute=$(pwd)
            path_relative=$(echo $path_absolute | sed -e 's|'$base_path'|.|')
            git_branch=$(git rev-parse --abbrev-ref HEAD)
            git_commit=$(git rev-parse HEAD)
            remote_url=$(git config --get remote.origin.url)
            guess_path=${remote_url/:/\/}
        cd - > /dev/null

        # Checking if repo was installed manually
        # A repo is installed manually if it has not .git folder (it is not a repo)
        # or it has git repo but this repo is not present in grasp-manager repos folder
        if [ -d "${installed_repo_paths[$i]}/.git" -a -d "$repos_path/$guess_path" ] 
        then
            echo "$log_git_token ${installed_repo_types[$i]} ${installed_repo_names[$i]} $path_relative $guess_path $git_branch $git_commit" >> $logfile
        else
            if [ -d "${installed_repo_paths[$i]}/.git" ]
            then
                echo "Manually installed and not restorable: ${installed_repo_types[$i]}-${installed_repo_names[$i]} $path_relative $git_branch $git_commit" >> $logfile
            else
                echo "Manually installed and not restorable: ${installed_repo_types[$i]}-${installed_repo_names[$i]} $path_relative" >> $logfile
            fi
        fi
    done
}

# This counter is a global variable which allows to know how many times
# prepare_environment function is called. It is because only in first call
# write_log function is called
prepare_environment_counter="0"

# This function read a yaml configuration file and prepare the environment
# $1 environment name
# $2 optional. Default value: grasp-manager.yml . You can specify other file to be parsed
function prepare_environment () {
    compilation_environment=$1
    let prepare_environment_counter=$prepare_environment_counter+1
    counter="0"

    # Load and prepare settings
    file_to_parse="grasp-manager.yml"
    if [ "$2" != "" ]
    then
        file_to_parse=$2
    fi

    eval $(parse_yaml $file_to_parse conf_)

    # Load and prepare compilation flags
    #if [ "$compilation_environment" = "~" ]
    #then
    #    compilation_environment="default"
    #fi

    using_last_environment="no"
    if [ "$compilation_environment" = "" ] # if there is not environment defined we use the last used for compilation and we don't update extensions
    then
        # I set it to special reserved word
        compilation_environment=$last_environment
        using_last_environment="yes"
    fi

    # Checking if compilation_environment is defined like default
    # In this case there is no valid settings (it will be loaded after)
    # but the compilation_environment is defined in the file
    # here we check partial names
    flag_name=conf_${compilation_environment}
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_compilation
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi    
    flag_name=conf_${compilation_environment}_extensions
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_input
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_input_drivers
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_input_transformers
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_output
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_output_segment_functions
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_output_current_functions
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    flag_name=conf_${compilation_environment}_output_tile_functions
    flag_value=${!flag_name}
    if [ "$flag_value" != "" ]
    then
        let counter=$counter+1
    fi 
    
    # Loading all settings
    flag_name=conf_${compilation_environment}_compilation_compile_dependencies
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_compile_dependencies=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_use_sudo
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_use_sudo=$flag_value
        let counter=$counter+1
    fi

    flag_name=conf_${compilation_environment}_compilation_constants_set
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_constants_set=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_build
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_build=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_mpi
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_mpi=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_f90
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_f90=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_prefix
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_prefix=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_build_dir
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_build_dir=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_cc
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_cc=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_cxx
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_cxx=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_compilation_cmake_flags
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        compilation_cmake_flags=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_modules_dependencies
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        modules_dependencies=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_modules_gpgpu
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        modules_gpgpu=$flag_value
        enable_retrieval_module gpgpu
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_modules_kernels_modifier
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        modules_kernels_modifier=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_modules_internal_examples
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        modules_internal_examples=$flag_value
        let counter=$counter+1
    fi
    flag_name=conf_${compilation_environment}_modules_models
    flag_value=${!flag_name}
    if [ "$flag_value" != "" -a "$flag_value" != "~" ]
    then
        modules_models=$flag_value
        let counter=$counter+1
    fi

    compilation_flags="CONSTANTS_SET=$compilation_constants_set BUILD=$compilation_build MPI=$compilation_mpi F90=$compilation_f90 PREFIX=$compilation_prefix BUILD_DIR=$compilation_build_dir CC=$compilation_cc CXX=$compilation_cxx CMAKE_FLAGS=\"$compilation_cmake_flags\" ENABLE_RETRIEVAL_MODULES=$retrieval_modules_enabled"

    env_split_by_underscores="${compilation_environment//[^_]}"
    env_number_of_underscores="${#env_split_by_underscores}"

    # Extension read in configuration file
    extension_drivers=$(set | grep ^conf_${compilation_environment}_extensions_input_drivers_ | cut -f1 -d= | cut -d_ -f$((6+$env_number_of_underscores))-99)
    extension_transformers=$(set | grep ^conf_${compilation_environment}_extensions_input_transformers_ | cut -f1 -d= | cut -d_ -f$((6+$env_number_of_underscores))-99)
    extension_segment_functions=$(set | grep ^conf_${compilation_environment}_extensions_output_segment_functions_ | cut -f1 -d= | cut -d_ -f$((7+$env_number_of_underscores))-99)
    extension_current_functions=$(set | grep ^conf_${compilation_environment}_xtensions_output_current_functions_ | cut -f1 -d= | cut -d_ -f$((7+$env_number_of_underscores))-99)
    extension_tile_functions=$(set | grep ^conf_${compilation_environment}_extensions_output_tile_functions_ | cut -f1 -d= | cut -d_ -f$((7+$env_number_of_underscores))-99)
    requested_kernels=$(set | grep ^conf_${compilation_environment}_kernels_ | cut -f1 -d= | cut -d_ -f$((4+$env_number_of_underscores))-99)
    requested_constants_sets=$(set | grep ^conf_${compilation_environment}_constants_sets_ | cut -f1 -d= | cut -d_ -f$((5+$env_number_of_underscores))-99)

    # If there is not assignation with current compilation environment
    # and it is not default environment means that the environment does not exist
    if [ "$counter" = "0" -a \
         "$compilation_environment" != "" -a \
         "$compilation_environment" != "~" -a \
         "$extension_drivers" == "" -a \
         "$extension_transformers" == "" -a \
         "$extension_segment_functions" == "" -a \
         "$extension_current_functions" == "" -a \
         "$extension_tile_functions" == "" -a \
         "$requested_kernels" == "" -a \
         "$requested_constants_sets" == "" ]
    then
        if [ "$using_last_environment" = "no" ]
        then
            gmecho "ERROR: '$compilation_environment' environment does not exist"
            exit 0
        else
            gmecho "WARNING: Last compilation environment ('$compilation_environment') does not exist and default values will be used."
            compilation_environment="~"
        fi
    fi

    # Total extensions the user wants to have
    drivers=$core_drivers" "$extension_drivers
    transformers=$core_transformers" "$extension_transformers
    segment_functions=$core_segment_functions" "$extension_segment_functions
    current_functions=$core_current_functions" "$extension_current_functions
    tile_functions=$core_tile_functions" "$extension_tile_functions
    kernels=$core_kernels" "$requested_kernels
    constants_sets=$core_constants_sets" "$requested_constants_sets

    # Current list of available extensions 
    installed_drivers=$(ls -F src/input/drivers | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1)
    installed_transformers=$(ls -F src/input/transformers | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1)
    installed_segment_functions=$(ls -F src/output/segment_functions | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1)
    installed_current_functions=$(ls -F src/output/current_functions | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1)
    installed_tile_functions=$(ls -F src/output/tile_functions | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1)
    installed_kernels=$(ls -F src/retrieval/internal_files | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1 |  grep -v ABS_GAS_DATA)
    installed_constants_sets=$(ls -F src/retrieval/constants_set | grep -e/ -e@ | cut -d/ -f1 | cut -d@ -f1 )

    # Current list of available extensions which don't come from core code. Current list of installed extensions
    installed_extension_drivers=$(not_in "$installed_drivers" "$core_drivers")
    installed_extension_transformers=$(not_in "$installed_transformers" "$core_transformers")
    installed_extension_segment_functions=$(not_in "$installed_segment_functions" "$core_segment_functions")
    installed_extension_current_functions=$(not_in "$installed_current_functions" "$core_current_functions")
    installed_extension_tile_functions=$(not_in "$installed_tile_functions" "$core_tile_functions")
    installed_external_kernels=$(not_in "$installed_kernels" "$core_kernels")
    installed_external_constants_sets=$(not_in "$installed_constants_sets" "$core_constants_sets")

    installed_repo_names=()
    installed_repo_types=()
    installed_repo_paths=()
    for driver in $installed_extension_drivers
    do
        installed_repo_names+=($driver)
        installed_repo_types+=("driver")
        installed_repo_paths+=("src/input/drivers/$driver")        
    done
    for transformer in $installed_extension_transformers
    do
        installed_repo_names+=($transformer)
        installed_repo_types+=("transformer")
        installed_repo_paths+=("src/input/transformers/$transformer") 
    done
    for segment_function in $installed_extension_segment_functions
    do
        installed_repo_names+=($segment_function)
        installed_repo_types+=("segment_function")
        installed_repo_paths+=("src/output/segment_functions/$segment_function")
    done
    for current_function in $installed_extension_current_functions
    do
        installed_repo_names+=($current_function)
        installed_repo_types+=("current_function")
        installed_repo_paths+=("src/output/current_functions/$current_function")
    done
    for tile_function in $installed_extension_tile_functions
    do
        installed_repo_names+=($tile_function)
        installed_repo_types+=("tile_function")
        installed_repo_paths+=("src/output/tile_functions/$tile_function")
    done
    for kernel in $installed_external_kernels
    do
        installed_repo_names+=($kernel)
        installed_repo_types+=("kernels")
        installed_repo_paths+=("src/retrieval/internal_files/$kernel")
    done
    for constants_set in $installed_external_constants_sets
    do
        installed_repo_names+=($constants_set)
        installed_repo_types+=("constants_set")
        installed_repo_paths+=("src/retrieval/constants_set/$constants_set")
    done

    # Modules

    # dependencies
    if [ -d $base_path/dependencies ]
    then
        installed_repo_names+=("dependencies")
        installed_repo_types+=("module")
        installed_repo_paths+=("dependencies")
    fi
    if [ -d $base_path/src/retrieval/gpgpu ]
    then
        installed_repo_names+=("gpgpu")
        installed_repo_types+=("module")
        installed_repo_paths+=("src/retrieval/gpgpu")
    fi
    if [ -d $base_path/src/retrieval/kernels_modifier ]
    then
        installed_repo_names+=("kernels_modifier")
        installed_repo_types+=("module")
        installed_repo_paths+=("src/retrieval/kernels_modifier")
    fi
    if [ -d $base_path/internal_examples ]
    then
        installed_repo_names+=("internal_examples")
        installed_repo_types+=("module")
        installed_repo_paths+=("internal_examples")
    fi
    if [ -d $base_path/src/retrieval/models ]
    then
        installed_repo_names+=("models")
        installed_repo_types+=("module")
        installed_repo_paths+=("src/retrieval/models")
    fi


    # Related to install and uninstall repos

    # List of extensions that it is necessary to install (they are not installed but they are in configuration file) 
    drivers_to_install=$(not_in "$extension_drivers" "$installed_drivers")
    transformers_to_install=$(not_in "$extension_transformers" "$installed_transformers")
    segment_functions_to_install=$(not_in "$extension_segment_functions" "$installed_segment_functions")
    current_functions_to_install=$(not_in "$extension_current_functions" "$installed_current_functions")
    tile_functions_to_install=$(not_in "$extension_tile_functions" "$installed_tile_functions")
    kernels_to_install=$(not_in "$requested_kernels" "$installed_kernels")
    constants_sets_to_install=$(not_in "$requested_constants_sets" "$installed_constants_sets")

    # List of extensions that it is necessary to uninstall (they are installed but not in configuration)
    drivers_to_remove=$(not_in "$installed_drivers" "$drivers")
    transformers_to_remove=$(not_in "$installed_transformers" "$transformers")
    segment_functions_to_remove=$(not_in "$installed_segment_functions" "$segment_functions")
    current_functions_to_remove=$(not_in "$installed_current_functions" "$current_functions")
    tile_functions_to_remove=$(not_in "$installed_tile_functions" "$tile_functions")     
    kernels_to_remove=$(not_in "$installed_kernels" "$kernels")
    constants_sets_to_remove=$(not_in "$installed_external_constants_sets" "$constants_sets")

    repos_to_install_names=()
    repos_to_install_types=()
    repos_to_install_paths=()
    repos_to_install_urls=()

    if [ ! -d $base_path/dependencies -a "$modules_dependencies" = "true" ]
    then
        repos_to_install_names+=("dependencies")
        repos_to_install_types+=("module")
        repos_to_install_paths+=("dependencies") 
        repos_to_install_urls+=($default_git_server_module_dependencies)
    fi
    if [ ! -d $base_path/src/retrieval/gpgpu -a "$modules_gpgpu" = "true" ]
    then
        repos_to_install_names+=("gpgpu")
        repos_to_install_types+=("module")
        repos_to_install_paths+=("src/retrieval/gpgpu") 
        repos_to_install_urls+=($default_git_server_module_gpgpu)
    fi
    if [ ! -d $base_path/src/retrieval/kernels_modifier -a "$modules_kernels_modifier" = "true" ]
    then
        repos_to_install_names+=("kernels_modifier")
        repos_to_install_types+=("module")
        repos_to_install_paths+=("src/retrieval/kernels_modifier") 
        repos_to_install_urls+=($default_git_server_module_kernels_modifier)
    fi
    if [ ! -d $base_path/internal_examples -a "$modules_internal_examples" = "true" ]
    then
        repos_to_install_names+=("internal_examples")
        repos_to_install_types+=("module")
        repos_to_install_paths+=("internal_examples") 
        repos_to_install_urls+=($default_git_server_module_internal_examples)
    fi
    if [ ! -d $base_path/src/retrieval/models -a "$modules_models" = "true" ]
    then
        repos_to_install_names+=("models")
        repos_to_install_types+=("module")
        repos_to_install_paths+=("src/retrieval/models") 
        repos_to_install_urls+=($default_git_server_module_models)
    fi
    
    for kernel in $kernels_to_install
    do
        repos_to_install_names+=($kernel)
        repos_to_install_types+=("kernels")
        repos_to_install_paths+=("src/retrieval/internal_files/$kernel")       

        url_stored=conf_${compilation_environment}_kernels_$kernel
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_kernels_prefix$kernel.git
        fi
        repos_to_install_urls+=($url) 
    done
    for constants_set in $constants_sets_to_install
    do
        repos_to_install_names+=($constants_set)
        repos_to_install_types+=("constants_set")
        repos_to_install_paths+=("src/retrieval/constants_set/$constants_set")       

        url_stored=conf_${compilation_environment}_constants_sets_$constants_set
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_constants_set_prefix$constants_set.git
        fi
        repos_to_install_urls+=($url) 
    done
    for driver in $drivers_to_install
    do
        repos_to_install_names+=($driver)
        repos_to_install_types+=("driver")
        repos_to_install_paths+=("src/input/drivers/$driver")       

        url_stored=conf_${compilation_environment}_extensions_input_drivers_$driver
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_driver_prefix$driver.git
        fi
        repos_to_install_urls+=($url) 
    done
    for transformer in $transformers_to_install
    do
        repos_to_install_names+=($transformer)
        repos_to_install_types+=("transformer")
        repos_to_install_paths+=("src/input/transformers/$transformer") 

        url_stored=conf_${compilation_environment}_extensions_input_transformers_$transformer
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_transformer_prefix$transformer.git
        fi
        repos_to_install_urls+=($url) 
    done
    for segment_function in $segment_functions_to_install
    do
        repos_to_install_names+=($segment_function)
        repos_to_install_types+=("segment_function")
        repos_to_install_paths+=("src/output/segment_functions/$segment_function")

        url_stored=conf_${compilation_environment}_extensions_output_segment_functions_$segment_function
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_segment_function_prefix$segment_function.git
        fi
        repos_to_install_urls+=($url) 
    done
    for current_function in $current_functions_to_install
    do
        repos_to_install_names+=($current_function)
        repos_to_install_types+=("current_function")
        repos_to_install_paths+=("src/output/current_functions/$current_function")

        url_stored=conf_${compilation_environment}_extensions_output_current_functions_$current_function
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_current_function_prefix$current_function.git
        fi
        repos_to_install_urls+=($url) 
    done
    for tile_function in $tile_functions_to_install
    do
        repos_to_install_names+=($tile_function)
        repos_to_install_types+=("tile_function")
        repos_to_install_paths+=("src/output/tile_functions/$tile_function")

        url_stored=conf_${compilation_environment}_extensions_output_tile_functions_$tile_function
        url=${!url_stored}
        if [ "$url" = "~" ]
        then
             url=$default_git_server:$default_git_tile_function_prefix$tile_function.git
        fi
        repos_to_install_urls+=($url) 
    done

    repos_to_uninstall_names=()
    repos_to_uninstall_types=()
    repos_to_uninstall_paths=()

    if [ -d $base_path/dependencies -a "$modules_dependencies" = "false" ]
    then
        repos_to_uninstall_names+=("dependencies")
        repos_to_uninstall_types+=("module")
        repos_to_uninstall_paths+=("dependencies") 
    fi
    if [ -d $base_path/src/retrieval/gpgpu -a "$modules_gpgpu" = "false" ]
    then
        repos_to_uninstall_names+=("gpgpu")
        repos_to_uninstall_types+=("module")
        repos_to_uninstall_paths+=("src/retrieval/gpgpu") 
    fi
    if [ -d $base_path/src/retrieval/kernels_modifier -a "$modules_kernels_modifier" = "false" ]
    then
        repos_to_uninstall_names+=("kernels_modifier")
        repos_to_uninstall_types+=("module")
        repos_to_uninstall_paths+=("src/retrieval/kernels_modifier") 
    fi
    if [ -d $base_path/internal_examples -a "$modules_internal_examples" = "false" ]
    then
        repos_to_uninstall_names+=("internal_examples")
        repos_to_uninstall_types+=("module")
        repos_to_uninstall_paths+=("internal_examples") 
    fi
    if [ -d $base_path/src/retrieval/models -a "$modules_models" = "false" ]
    then
        repos_to_uninstall_names+=("models")
        repos_to_uninstall_types+=("module")
        repos_to_uninstall_paths+=("src/retrieval/models") 
    fi

    for kernel in $kernels_to_remove
    do
        repos_to_uninstall_names+=($kernel)
        repos_to_uninstall_types+=("kernels")
        repos_to_uninstall_paths+=("src/retrieval/internal_files/$kernel")       
    done
    for constants_set in $constants_sets_to_remove
    do
        repos_to_uninstall_names+=($constants_set)
        repos_to_uninstall_types+=("constants_set")
        repos_to_uninstall_paths+=("src/retrieval/constants_set/$constants_set")       
    done
    for driver in $drivers_to_remove
    do
        repos_to_uninstall_names+=($driver)
        repos_to_uninstall_types+=("driver")
        repos_to_uninstall_paths+=("src/input/drivers/$driver")       
    done
    for transformer in $transformers_to_remove
    do
        repos_to_uninstall_names+=($transformer)
        repos_to_uninstall_types+=("transformer")
        repos_to_uninstall_paths+=("src/input/transformers/$transformer") 
    done
    for segment_function in $segment_functions_to_remove
    do
        repos_to_uninstall_names+=($segment_function)
        repos_to_uninstall_types+=("segment_function")
        repos_to_uninstall_paths+=("src/output/segment_functions/$segment_function")
    done
    for current_function in $current_functions_to_remove
    do
        repos_to_uninstall_names+=($current_function)
        repos_to_uninstall_types+=("current_function")
        repos_to_uninstall_paths+=("src/output/current_functions/$current_function")
    done
    for tile_function in $tile_functions_to_remove
    do
        repos_to_uninstall_names+=($tile_function)
        repos_to_uninstall_types+=("tile_function")
        repos_to_uninstall_paths+=("src/output/tile_functions/$tile_function")
    done

                     

    if [ "$prepare_environment_counter" = "1" ]
    then
        write_log
        if [ "$using_last_environment" = "yes" ]
        then
            gmecho "ENVIRONMENT NOT SPECIFIED: "
            gmecho "   Using $compilation_environment environment for compilation configuration and current installed extensions as extensions"
        fi
    fi
}

function enable_retrieval_module () {
    if [ "$retrieval_modules_enabled" != "" ]
    then
        retrieval_modules_enabled+=","
    fi
    retrieval_modules_enabled+=$1
}

# This function based in what the users want (read from setting file) update
# available extension for the user and prepare again the environment 
function update_extensions () {
    uninstall_extensions
    install_extensions
    prepare_environment "$compilation_environment"
}

function git_apply () {
    i=0

    # obtaining backup folder with absolute path

    if [ "$?" != "0" -o ! -d "$1" ]
    then
        gmecho_title "ERROR: back up path does not exist"
    fi

    backup_path=$1

    # Removing all extensions.
    tmp_yaml_file=$home_manager"/tmp-grasp-manager.yml" # Preparing the temporal file
    echo -e "\n$last_environment: ~\n" > $tmp_yaml_file  # Printing in temporal file
    prepare_environment "$last_environment" "$tmp_yaml_file" # Call to install action
    rm $tmp_yaml_file # Removing temporal file
    uninstall_extensions # installing extension


    token_length=$(echo -ne $log_git_token | wc -c)
    cat -v "$backup_path/$logfilename" | while read line
    do
        let i=$i+1
        if [ "${line:0:$token_length}" = "$log_git_token" ]
        then
            repo_type=$(echo $line|cut -d ' ' -f2)
            repo_name=$(echo $line|cut -d ' ' -f3)
            target_folder=$(echo $line|cut -d ' ' -f4)
            origin_folder=$(echo $line|cut -d ' ' -f5)
            origin_folder=$repos_path/$origin_folder
            branch_name=$(echo $line|cut -d ' ' -f6)
            commit=$(echo $line|cut -d ' ' -f7)

            # Checking if the origin folder is available, otherwise the process is stopped
            if [ ! -d $origin_folder ]
            then
                # Removing all extensions.
                tmp_yaml_file=$home_manager"/tmp-grasp-manager.yml" # Preparing the temporal file
                echo -e "\n$last_environment: ~\n" > $tmp_yaml_file  # Printing in temporal file
                prepare_environment "$last_environment" "$tmp_yaml_file" # Call to install action
                rm $tmp_yaml_file # Removing temporal file
                uninstall_extensions # installing extension

                gmecho_title "An error has been detected trying to install $repo_type $repo_name. Process aborted".
                exit -1
            fi

            # stablish the code
            cp -a $origin_folder $target_folder

            # We stablish the commit
            # We set the branch name if it exist
            # if once we stablish branch name is not same commit, we move to the commit
            cd $target_folder

            if [ "$branch_name" = "HEAD" ]
            then
                git checkout $commit
            else
                git checkout $branch_name
                current_commit=$(git rev-parse HEAD)
                if [ "$current_commit" != "$commit" ]
                then
                    git checkout $commit
                fi
            fi
            
            # Finally we apply the patch if it exist
            if [ -f $backup_path/$repo_type-$repo_name.patch ]
            then
                git apply $backup_path/$repo_type-$repo_name.patch
            fi
            cd - > /dev/null
        fi
    done 

    # Applying changes to core
    core_line=$(cat $backup_path/$logfilename | grep "^#core_restore_info")
    
    branch_name=$(echo $line|cut -d ' ' -f3)
    commit=$(echo $line|cut -d ' ' -f4)

    if [ "$branch_name" = "HEAD" ]
    then
        git checkout $commit
    else
        git checkout $branch_name
        current_commit=$(git rev-parse HEAD)
        if [ "$current_commit" != "$commit" ]
        then
            git checkout $commit
        fi
    fi

    if [ -f $backup_path/core.patch ]
    then
        git apply $backup_path/core.patch
    fi
    
}


function create_extension_function () {
    extension_type=$1
    extension_name=$2

    if [ "$3" != "" ]
    then
        extension_url=$3
    else
        if [ "$extension_type" = "driver" ]
        then
            extension_url=$default_git_server:$default_git_driver_prefix$extension_name
        fi
        if [ "$extension_type" = "transformer" ]
        then
            extension_url=$default_git_server:$default_git_transformer_prefix$extension_name
        fi
        if [ "$extension_type" = "segment_function" ]
        then
            extension_url=$default_git_server:$default_git_segment_function_prefix$extension_name
        fi
        if [ "$extension_type" = "current_function" ]
        then
            extension_url=$default_git_server:$default_git_current_function_prefix$extension_name
        fi
        if [ "$extension_type" = "tile_function" ]
        then
            extension_url=$default_git_server:$default_git_tile_function_prefix$extension_name
        fi        
    fi

    # Setting bsae path
    base_extension_path=""
    io=""
    main_function_name=""
    main_function_interface=""
    include_files_c=""
    include_files_h=""
    if [ "$extension_type" = "driver" ]
    then
        main_function_structure_field="get_segment"
    else
        main_function_structure_field="function"
    fi

    if [ "$extension_type" = "driver" ]
    then
        base_extension_path="src/input/drivers"
        io="input"
        main_function_name="get_segment"
        main_function_interface="grasp_settings *settings, grasp_segment_t *segment, int col, int row, int itime"
        include_files_c='#include "../../../input/grasp_input.h"\n#include "../../../input/grasp_input_functions.h"'
        include_files_h='#include "../../../settings/grasp_settings.h"\n#include "../../../input/grasp_input.h"\n#include "../../../input/grasp_input_functions.h"'
    fi
    if [ "$extension_type" = "transformer" ]
    then
        base_extension_path="src/input/transformers"
        io="input"
        main_function_name="function"
        main_function_interface="grasp_settings *settings,grasp_segment_t *segment"
        include_files_c='#include "../../../input/grasp_input.h"\n#include "../../../input/grasp_input_functions.h"'
        include_files_h='#include "../../../settings/grasp_settings.h"\n#include "../../../input/grasp_input.h"\n#include "../../../input/grasp_input_functions.h"'
    fi
    if [ "$extension_type" = "segment_function" ]
    then
        base_extension_path="src/output/segment_functions"
        io="output"
        main_function_name="process"
        main_function_interface="grasp_output_stream *stream, grasp_settings *settings,grasp_segment_t *segment,output_segment_general *output,grasp_tile_description_t *tile_description,int icol,int irow,int itime"
        include_files_h='#include "../../grasp_output_stream_t.h"\n#include "../../../settings/grasp_settings.h"\n#include "../../../input/grasp_input.h"\n#include "../../grasp_output.h"'
    fi
    if [ "$extension_type" = "current_function" ]
    then
        base_extension_path="src/output/current_functions"
        io="output"
        main_function_name="process"
        main_function_interface="grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *current, output_segment_general *output, grasp_results_t *results, grasp_tile_description_t *tile_description,int icol,int irow,int itime"
        include_files_h='#include "../../grasp_output_stream_t.h"\n#include "../../../settings/grasp_settings.h"\n#include "../../../input/grasp_input.h"\n#include "../../grasp_output.h"'
    fi
    if [ "$extension_type" = "tile_function" ]
    then
        base_extension_path="src/output/tile_functions"
        io="output"
        main_function_name="process"
        main_function_interface="grasp_output_stream *stream, grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results"
        include_files_h='#include "../../grasp_output_stream_t.h"\n#include "../../../settings/grasp_settings.h"\n#include "../../../input/grasp_input.h"\n#include "../../grasp_output.h"'
    fi   

    extension_path=$base_extension_path"/"$extension_name
    mkdir -p $extension_path

    # Prepare readme file
    echo -ne "# $extension_name $extension_type\n## GRASP EXTENSION\n\n\nThis README file has to be completed" > $extension_path/README.md

    # Initial code files
    # Cmake file
    echo -ne "add_library(${io}_${extension_type}_${extension_name} grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}.c)
target_link_libraries(${io}_${extension_type}_${extension_name} grasp_settings)
" > $extension_path/CMakeLists.txt
# this line has been removed since it has been moved to input module (cmake file inside extensions)
# set_property(GLOBAL APPEND PROPERTY $(tr '[a-z]' '[A-Z]' <<< $io)_$(tr '[a-z]' '[A-Z]' <<< $extension_type)S ${extension_name})


    # Code files
echo -ne "
#include <stdio.h>
#include \"grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}.h\"
${include_files_c}

grasp_${io}_${extension_type}_t grasp_${io}_${extension_type}_${extension_name}(){
    grasp_${io}_${extension_type}_t x;
    
    x.init=grasp_${io}_${extension_type}_${extension_name}_init;
    x.${main_function_structure_field}=grasp_${io}_${extension_type}_${extension_name}_${main_function_name};
    x.close=grasp_${io}_${extension_type}_${extension_name}_close;

    return x;
}

int grasp_${io}_${extension_type}_${extension_name}_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_${io}_${extension_type}_${extension_name}_close(void){
    return 0;
}

int grasp_${io}_${extension_type}_${extension_name}_${main_function_name}(${main_function_interface}){
    return 0;
}

grasp_settings_parameter_array *grasp_${io}_${extension_type}_settings_${extension_name}(grasp_settings *settings){
    return NULL;
}
" > $extension_path/grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}.c

echo -ne "
#ifndef $(tr '[a-z]' '[A-Z]' <<< grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}_h)
#define	$(tr '[a-z]' '[A-Z]' <<< grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}_h)

#ifdef	__cplusplus
extern "C" {
#endif

${include_files_h}   

grasp_${io}_${extension_type}_t grasp_${io}_${extension_type}_${extension_name}();

int grasp_${io}_${extension_type}_${extension_name}_init(grasp_settings *settings, grasp_tile_description_t *input_information);

int grasp_${io}_${extension_type}_${extension_name}_close(void);

int grasp_${io}_${extension_type}_${extension_name}_${main_function_name}(${main_function_interface});

grasp_settings_parameter_array *grasp_${io}_${extension_type}_settings_${extension_name}(grasp_settings *settings);

#ifdef	__cplusplus
}
#endif

#endif	/* $(tr '[a-z]' '[A-Z]' <<< grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}_h) */
" > $extension_path/grasp_${io}_$( echo $extension_type | cut -f1 -d_)_${extension_name}.h

    cd $extension_path
    
    # Init git repository
    git init
    git add -A .
    git commit -m "Initial commit"
    
    # Configure branches
    git checkout -b dev

    # Synch with server
    gmecho "If you want to pull the code to a server, please create first the repository $extension_url"
   
    while : ; do # Looping to upload the code until it works
        gmecho "Do you want to continue pulling the code [Y/n]:"
        read answer
        git_response="0"
        if [ "$answer" = "" -o "$answer" = "Y" -o "$answer" = "y" ] 
        then
            git remote add origin $extension_url
            git push origin master
            let git_response=$git_response+$?
            git push origin dev
            let git_response=$git_response+$?
        else
            git_response="0"
        fi
        if [ "$git_response" = "0" ]
        then
            break;
        else
            gmecho "There was an error uploading the code. Have you created the repository in git server?"
        fi
    done

    cd - > /dev/null
}

function parse_yaml_error {
    echo -e $1
    exit -1
}

function parse_yaml {
   # Base code Stefan Farestam: http://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_-]*' fs=$(echo @|tr @ '\034')
   tr -d $'\r'< $1 | 
   sed -ne "s|^\($s\):|\1|" \
        -e 's/#.*$//' \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"   |
   awk -F$fs '{
      indent = length($1)/4;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      gsub(/[ \t]+$/, "", $3)
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
      if ($2 ~ '/-/') { printf("parse_yaml_error \"Error parsing YAML file: It constains an dash (-) symbol in left side. \
                                 \\nIt is not allowed in current version of grasp-manager. \
                                 \\nIf you have a definition like: \
                                 \\nenv_name: \
                                 \\n     parameter-name: value \
                                 \\nReplace parameter-name by parameter_name. \
                                 \\nSpecifically the error was found in %s\"", $2);exit -1; }
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

function is_in () {
    [[ $2 =~ $1 ]] && echo "1" || echo "0"
}

function not_in () {
    for i in $1
    do
        found="0"
        for j in $2
        do
            if [ "$i" = "$j" ]
            then
                found="1"
                break
            fi
        done
        if [ "$found" = "0" ]
        then
            echo $i
        fi
    done
}

# This function select a extension version given core version
# If extension does not have version numbers it will return master
# if extension has not exactly the version number it will return the biggest version smaller than given 
function select_extension_version () {

    search_version=$1
    selected_version=""
    # First we check if search_version is a branch o tag
    # if is a tag is_version is 'yes'
    # to know it, we go to core folder (cd ..) and we check if
    # current version is in the list of tags.
    cd .. > /dev/null
    commit_type="branch_name"
    git tag | grep $search_version > /dev/null
    if [ "$?" = "0" ]
    then
        commit_type="version_number"
    fi
    # We need to know if the search version is a commit name
    # if it is a commit name we have to move to best version or to master
    git rev-list --all --remotes | grep ^$search_version > /dev/null
    if [ "$?" = "0" ]
    then
        commit_type="commit"
        # search_version will be an equivalent to the biggest version which contain the commit
        # It is like say that this commit is compatible with this version
        search_version=$(git tag --contains $search_version | grep ^v | sort -t. -k1.2,1nr -k2,2nr -k3,3nr | head -n 1)
    fi
    cd - > /dev/null
    
    extension_list_branches=$((git branch -r --no-color | cut -c3-999 | cut -d/ -f2-99 | grep -v ^HEAD && git branch --no-color | cut -c3-99) | sort -u)

    if [ "$commit_type" = "branch_name" ] 
    then
        # if it is a name of branch we move the extension repository to 
        # name of branch, if it does not exist to dev, and finally, if it does not
        # exist we move it to master
        
        is_branch_available=$(is_in "$search_version" "$extension_list_branches")
        
        if [ "$is_branch_available" = "1" ] # Branch name exists, we return it
        then
            selected_version=$search_version
        else
            is_branch_available=$(is_in "dev" "$extension_list_branches")
            if [ "$is_branch_available" = "1" ] # If we branch name does not exist we try we dev
            then
                selected_version="dev"
            else # we return master
                select_version="master"
            fi
        fi
    else 
        #  "$commit_type" = "version_number" or "$commit_type" = "commit" 

        version_prepared=$search_version".99"
        # we prepare a ordered list of versions where after $version_prepared is our target version
        versions=$((echo $version_prepared &&  git tag) | grep ^v | sort -t. -k1.2,1nr -k2,2nr -k3,3nr -k4,4r | grep ^$search_version)
        found=0

        for version in $versions
        do
            if [ $found = 1 ] # We have found next version of our lure (because it was in [ "$version" = "$version_prepared" ] so it is our target version
            then
                selected_version=$version
                break
            fi
            if [ "$version" = "$version_prepared" ] # We found our lure so next version is the good one
            then
                found=1
            fi
        done
        if [ "$selected_version" = "" ] # If we have not found version number number valid we return branches. Dev and if does not exist, master
        then
            is_branch_available=$(is_in "dev" "$extension_list_branches")
            if [ "$is_branch_available" = "1" ] # If we branch name does not exist we try we dev
            then
                selected_version="dev"
            else # we return master
                select_version="master"
            fi
        fi
    fi
    
    echo $selected_version
}

function register_password {
    gmecho_title "The action you want to execute require sudo access"
   
    while : ; do  # Ask the password until it is valid
        gmecho_title "Please, if you want to continue provide the password:"
        read -s password
        echo $password | sudo -S ls > /dev/null 2> /dev/null
        if [ $? = 0 ] 
        then
            break
        fi
    done

    # Set sudo variable
    sudo="echo $password | sudo -S "
}

##########################################
####  Definition of action functions  ####
##########################################

function update-grasp () {
    prepare_environment "$1"
    if [ "$compilation_use_sudo" = "true" ]
    then
        register_password
    fi
    backup
    update_extensions
    git_fetch
    git_checkout "$last_version_tag"
    git_update_from_remote
    make_dependencies
    make_compile
    make_install
}

function update-grasp-to () {
    prepare_environment "$2"
    if [ "$compilation_use_sudo" = "true" ]
    then
        register_password
    fi
    backup
    update_extensions
    git_fetch
    git_checkout "$1"
    git_update_from_remote
    make_dependencies
    make_compile
    make_install
}

function update-grasp-to-dev () {
    prepare_environment "$1"
    if [ "$compilation_use_sudo" = "true" ]
    then
        register_password
    fi
    backup
    update_extensions
    git_fetch
    git_checkout "dev"
    git_update_from_remote
    make_dependencies
    make_compile
    make_install
}

function make () {
    prepare_environment "$1"
    make_compile
}

function clean () {
    make_clean
}

function purge () {
    if [ "$1" != "all-backups" -a "$1" != "unused-repos" -a "$1" != "all" -a "$1" != "manual-installed"  -a "$1" != "all-repos" -a "$1" != "all-installed-repos" ]
    then
        echo "ERROR: You need to specify what do you want to purge. Valid options are all, all-backups and all"
        exit 0
    fi

    manager_purge "$1"
}

function pull () {
    prepare_environment "$1"
    backup
    update_extensions
    git_fetch
    git_checkout "$current_version"
    git_update_from_remote    
}

function checkout () {
    if [ "$1" = "" ]
    then
        echo "ERROR: commit, tag or branch name is mandatory. Please run ./grasp-manager.sh checkout commit"
        exit 0
    fi
    # Checking if branch exists
    verified="no"
    git rev-parse --verify $1 > /dev/null 2>&1
    if [ "$?" = "0" ]
    then
        verified="yes"
    else
        git rev-parse --verify origin/$1 > /dev/null 2>&1
        if [ "$?" = "0" ]
        then
            verified="yes"
        fi
    fi

    if [ "$verified" = "yes" ]
    then
        prepare_environment "$2"
        backup
        update_extensions
        git_fetch
        git_checkout "$1"
        git_update_from_remote
    else
        gmecho_title "ERROR: $1 is not a valid git reference in current GRAPS repository"
    fi
}

function apply () {
    if [ "$1" = "" ]
    then
        echo "ERROR: backup-path is mandatory. Please run ./grasp-manager.sh apply $bk_folder_path[backup-folder]"
        exit 0
    fi
    #prepare_environment
    backup
    git_apply "$bk_folder_path$1"        
}

function rollback () {
    #prepare_environment
    backup
    # Obtaining backup path
    git_apply "$bk_folder_path$last_backup"      
}

function commit () {
    if [ "$1" = "" ]
    then
        echo "ERROR: message is mandatory. Please run ./grasp-manager.sh commit messages"
        exit 0
    fi
    if [ "$1" = "-m" ]
    then
        echo "ERROR: Remember that grasp-manager is not git. -m option is not necessary, just first argument is the message"
        exit 0
    fi
    prepare_environment
    git_add
    git_commit "$1"
    #git_fetch
    #git_update_from_remote
    #git_push
}

function status () {
    prepare_environment 
    #git_fetch
    git_status
}

function tag () {
    if [ "$1" = "" ]
    then
        echo "ERROR: tag number is mandatory. Please run ./grasp-manager.sh tag tag_number environment"
        exit 0
    fi
    prepare_environment
    backup
    #git_checkout "master"
    git_tag "$1"
}

function merge () {
    if [ "$1" = "" ]
    then
        echo "ERROR: branch name is mandatory. Please run ./grasp-manager.sh merge 'branch_name'"
        exit 0
    fi
    prepare_environment
    git_merge "$1"
}

function install_extension () {
    if [ "$1" != "driver" -a "$1" != "transformer" -a "$1" != "segment_function" -a "$1" != "current_function" -a "$1" != "tile_function" -a "$1" != "module" -a "$1" != "kernels" -a "$1" != "constants_set" ]
    then
        echo "ERROR: Extension type has to be 'driver' or 'transformer' or 'segment_function' or 'current_function or 'tile_function' or 'module' or 'kernels' or 'constants_set'"
        exit 0
    fi

    manager_install_extension "$1" "$2" "$3"
}

function uninstall_extension () {
    if [ "$1" != "driver" -a "$1" != "transformer" -a "$1" != "segment_function" -a "$1" != "current_function" -a "$1" != "tile_function" -a "$1" != "module" -a "$1" != "kernels" -a "$1" != "constants_set" ]
    then
        echo "ERROR: Extension type has to be 'driver' or 'transformer' or 'segment_function' or 'current_function or 'tile_function' or 'module' or 'kernels' or 'constants_set'"
        exit 0
    fi

    manager_uninstall_extension "$1" "$2"
}

function branch () {
    if [ "$1" = "" ]
    then
        echo "ERROR: branch name is mandatory. Please run ./grasp-manager.sh branch 'branch_name'"
        exit 0
    fi
    prepare_environment
    git_branch "$1"
}

function push () {
    prepare_environment
    git_push
}

function create_extension () {
    if [ "$1" != "driver" -a "$1" != "transformer" -a "$1" != "segment_function" -a "$1" != "current_function" -a "$1" != "tile_function" ]
    then
        echo "ERROR: Extension type has to be 'driver' or 'transformer' or 'segment_function' or 'current_function or 'tile_function'"
        exit 0
    fi

    if [ "$2" = "" ]
    then
        echo "ERROR: argument two has to be the name of the extension"
        exit 0
    fi

    prepare_environment

    create_extension_function "$1" "$2" "$3"
}


#########################
##### Main function #####
#########################

if [ $# -eq 0 ]
then    
    help
    exit 0
fi

action=$1
shift

case "$action" in
    "update-grasp")
        update-grasp "$1"
        ;;
    "update-grasp-to")
        update-grasp-to "$1" "$2"
        ;;
    "update-grasp-to-dev")
        update-grasp-to-dev "$1"
        ;;
    "commit")
        commit "$1"
        ;;
    "tag")
        tag "$1"
        ;;
    "checkout")
        checkout "$1" "$2"
        ;;
    "status")
        status "$1"
        ;;
    "apply")
        apply "$1"
        ;;
    "rollback")
        rollback "$1"
        ;;
    "pull")
        pull "$1" 
        ;;
    "push")
        push "$1" 
        ;;
    "merge")
        merge "$1" 
        ;;
    "branch")
        branch "$1" 
        ;;
    "create-extension")
        create_extension "$1" "$2" "$3"
        ;;
    "make")
        make "$1"
        ;;
    "install-extension")
        install_extension "$1" "$2" "$3"
        ;;
    "uninstall-extension")
        uninstall_extension "$1" "$2"
        ;;
    "clean")
        clean
        ;;
    "purge")
        purge "$1"
        ;;
    *)
        echo "ERROR: Action not defined"
        help
        ;;
esac


