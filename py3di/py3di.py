import os
import sys
import logging
import argparse
import subprocess
from corewrapper.m_simulation import Simulation
from corewrapper.m_model import Modelconfiguration
from corewrapper.m_paths import Projectpaths


class Py3di(object):

    def __init__(self, ini_file_path, filepaths=None):

        project_root = ini_file_path.split('.')[0]
        project = ProjectManager()
        project.create_basic_layout(base_path=project_root)

        if sys.platform == 'win32':
            model_root = ini_file_path.split('\\')[0] + '/'
        else:
            model_root = ini_file_path.split('/')[0] + '/'
                    
        if not filepaths:
            project_root = project_root + '/' 
            input_path = project_root + 'input_generated/'
            log_path = project_root + 'log/'
            result_path = project_root + 'results/'
            admin_path = project_root + 'preprocessed/'
            rain_path = project_root + 'scenario/rainfall/'
            boundary_path = project_root + 'scenario/boundary_conditions/'
            external_input_path = project_root + 'data/'
            state_path = project_root + 'states/'
        else:
            'Blah' 

        project_paths = Projectpaths(model_root, project_root, input_path, log_path, 
                                    result_path, admin_path, rain_path, 
                                    boundary_path, external_input_path, 
                                    state_path)

        self.model = Modelconfiguration(ini_file_path, project_paths)


    def initiate_simulation(self):
        
        self.sim = Simulation(self.model)


class ProjectManager(object):

    def __init__(self):

        self.tree = None
        self.path_map = dict()

    def _populate_tree(self, tree_dict=None, base_path=None, create_dirs=True):
        """
        Recursively creates all directories in the tree.
        :param tree_dict: dict of dict (max depth 1000)
                          of directories to create
        :param base_path: starting path
        :return:
        """

        if tree_dict is None:
            tree_dict = self.tree

        for dir_name, subdir_name in tree_dict.items():
            if hasattr(subdir_name, 'items'):
                logger.debug("[DB] has subdirs, digging deeper...")
                current_path = os.path.join(base_path, dir_name)
                logger.debug("[DB] current_path is: {}".format(current_path))
                self._populate_tree(
                    subdir_name, base_path=current_path,
                    create_dirs=create_dirs)
            else:
                dest = os.path.join(base_path, dir_name)
                self.path_map[dir_name] = dest
                if create_dirs is False:
                    continue

                if not os.path.exists(dest):
                    try:
                        os.makedirs(dest)
                        logger.info(
                            "[+] Successfully created directory: "
                            "{dirname} ".format(dirname=dest)
                        )
                    except OSError:
                        logger.exception(
                            "[-] Could not create directory "
                            "{dir_name}!".format(dir_name=dest)
                        )
                else:
                    logger.info(
                        "[*] Directory already exists. Skipping creation..."
                    )

    def create_basic_layout(self, base_path=None, tree=None, create_dirs=True):
        """
        Wrapper around the ``_populate_tree`` method. Will create all
        directories as specified by the ``tree`` parameter.
        :param tree: dict of dict like so:
            {<dir_name_1>: {<subdir_1_name> : None}, {<subdir_2_name> : None},
            {<subdir_n_name> : {subsubdir_name_1: None}}, <dir_name_2>: None,
            ...}
        """
        self.tree = tree
        if self.tree is None:
            self.tree = {
                'modules': None,
                'input_generated': None,
                'scenario': {
                    'boundary_conditions': None,
                    'rainfall': None,
                    'initial_conditions': None,
                    'wind': None,
                    'laterals': None,
                    'controls': None,
                    },
                'log': None,
                'data': None,
                'results': None,
                'preprocessed': None,
                'states': None,
            }

        self._populate_tree(
            self.tree, base_path=base_path, create_dirs=create_dirs)
