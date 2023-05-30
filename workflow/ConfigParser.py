
from typing import List, Tuple, Optional
from pathlib import Path
import itertools
import pandas as pd

# How to add new parameters:
# 1. Add the parameter to DEFAULT_PARAMETERS
# 2. Add the parameter to DEFAULT_PARAMETERS_TYPES
# 3. If the parameter is a pseudo parameter (i.e. it's not used directly but converted into other parameters):
#    3.1 Add the parameter to PSEUDO_PARAM_TO_PARAM 
#    3.2 Add the conversion in the ConfigParser class at `self._dataset_pseudo_parameter_conversion` or
#        `self._selection_pseudo_parameter_conversion`
# 4. Adjust the scripts to support the new parameter (process_data.py, selection.py, etc.)

DEFAULT_PARAMETERS = {
    "dataset": {
        "processing": None, # options: [None, "lognorm"]
        "ct_key": "celltype", # only relevant if e.g. n_cts is specified
        "n_cts" : None,
        "cells_per_ct_n_seeds" : 1, # Pseudo parameter
        "cells_per_ct_seed": 0,
        "cells_per_ct": None,
    },
    "selection": {
        "n": 100,
        "ct_key": "celltype",
        "gene_key": None,
        "method_specific_processing": False, # wether to use the method specific normalization/filtering etc. 
    },
}

"""Some parameters are not used directly and instead converted into (multiple) values of other parameters.

For example, `cells_per_ct_n_seeds = 3` is converted into three configurations with `cells_per_ct_seed = 0`, 
`cells_per_ct_seed = 1`, and `cells_per_ct_seed = 2`.

"""
PSEUDO_PARAM_TO_PARAM = {
    "dataset": {
        "cells_per_ct_n_seeds": "cells_per_ct_seed",
    },
    "selection": {},
}
PARAM_TO_PSEUDO_PARAM = {
    "dataset" : {v: k for k,v in PSEUDO_PARAM_TO_PARAM["dataset"].items()},
    "selection" : {v: k for k,v in PSEUDO_PARAM_TO_PARAM["selection"].items()},
}

# Convert None to str (reading yamls interprets None as str and when saving the dataframes we want to keep None instead of empty fields)
DEFAULT_PARAMETERS = {
    "dataset": {k: str(v) if v is None else v for k,v in DEFAULT_PARAMETERS["dataset"].items()},
    "selection": {k: str(v) if v is None else v for k,v in DEFAULT_PARAMETERS["selection"].items()},
}

DEFAULT_PARAMETERS_TYPES = {
    "dataset": {
        "dataset": str, # this one is special as it's separate from the dataset_param
        "processing": str,
        "ct_key": str, # only relevant if e.g. n_cts is specified
        "n_cts" : int,
        "cells_per_ct_n_seeds" : int,
        "cells_per_ct_seed": int,
        "cells_per_ct": int,
    },
    "selection": {
        "method": str, # this one is special as it's separate from the selection_param 
        "n": int,
        "ct_key": str,
        "gene_key": str,
        "method_specific_processing": bool,
    },
}

DEFAULT_PARAMETERS_NON_PSEUDO = {
    "dataset": {k: v for k,v in DEFAULT_PARAMETERS["dataset"].items() if k not in PSEUDO_PARAM_TO_PARAM["dataset"]},
    "selection": {k: v for k,v in DEFAULT_PARAMETERS["selection"].items() if k not in PSEUDO_PARAM_TO_PARAM["selection"]},
}



class ConfigParser():
    
    def __init__(self, config: dict) -> None:
        """
        Parse the config file and generate the final file names
        
        We create three tables:
        1. Table of dataset configurations and their respective ids (1 row per unique configuration)
        2. Table of selection configurations and their respective selection and dataset ids (1 row per unique configuration)
        3. Overview table of all selections defined in each batch
           (1 row per selection, different batches can include the same selection)
        
        
        
        config: dict
            Dict of config yaml
        """
        
        ## Read yaml config file
        #with open(config, 'r') as stream:
        #    self.config = yaml.safe_load(stream) 
        #self.config = config
        
        # Convert pseudo parameters to parameters
        self.config = self._convert_pseudo_params_to_param(config.copy())
        
        # Set paths
        self.DATA_DIR = self.config['DATA_DIR']
        self.DATA_DIR_TMP = self.config['DATA_DIR_TMP']
        self.RESULTS_DIR = self.config['RESULTS_DIR']
        self.SAVE_METHOD_SPECIFIC_OUTPUT = self.config['SAVE_METHOD_SPECIFIC_OUTPUT']
        
        # Make dirs
        Path(self.DATA_DIR_TMP).mkdir(exist_ok=True)
        Path(self.RESULTS_DIR).mkdir(exist_ok=True)
        
        # File names of configuration tables
        self.data_params_file = Path(self.RESULTS_DIR, "data_parameters.csv")
        self.selection_params_file = Path(self.RESULTS_DIR, "selection_parameters.csv")
        self.selection_overview_file = Path(self.RESULTS_DIR, "selection_overview.csv")
        
        ## Reduce default hyperparameters to those that occur in the config
        ## TODO (not needed anymore with current solution): If configuration files of previous runs exist, check if there are additional hparams that were used previously
        #default_hparams = self._get_hyperparams_that_occur_in_config()
        # We just use all default hyperparameters since it's not that many
        default_hparams = DEFAULT_PARAMETERS_NON_PSEUDO
        
        # Dataset configurations
        dataset_param_combs = {}
        for batch, batch_dict in self.config['selections'].items():            
            dataset_param_combs[batch] = self._get_combinations_of_params(
                batch_dict, default_hparams,
                main_key = "datasets",
                main_key_param_name = "dataset",
                param_key = "dataset_param",
                default_param_key = "dataset",
            )            
            
        # Selection configurations
        selection_param_combs = {}
        for batch, batch_dict in self.config['selections'].items():
            selection_param_combs[batch] = self._get_combinations_of_params(
                batch_dict, default_hparams,
                main_key = "methods",
                main_key_param_name = "method",
                param_key = "selection_param",
                default_param_key = "selection",
            )
            
        # Load configuration files of previous to conserve old ids and add new ids respectively
        data_ids_to_config_old = self._load_config_table(self.data_params_file, param_group="dataset") if self.data_params_file.exists() else None
        selection_ids_to_config_old = self._load_config_table(self.selection_params_file, param_group="selection") if self.selection_params_file.exists() else None
        
        # Set dataset ids and add ids to dataset configurations
        data_ids_to_config, dataset_param_combs = self._get_ids_and_configs(
            dataset_param_combs, id_str="data_id", ids_to_config_old=data_ids_to_config_old
        )
        # Set selection ids and add ids to selection configurations
        selection_ids_to_config, selection_param_combs = self._get_ids_and_configs(
            selection_param_combs, id_str="selection_id", ids_to_config_old=selection_ids_to_config_old
        )
        
        self.dfs = {}
        
        # Data configurations table
        df_data = pd.DataFrame(data_ids_to_config.values())
        df_data["data_id"] = data_ids_to_config.keys()
        df_data = df_data.astype("object")
        self.dfs["data"] = df_data.set_index("data_id")
        
        # Selection configurations table
        df_selection = pd.DataFrame(selection_ids_to_config.values())
        df_selection["selection_id"] = selection_ids_to_config.keys()
        df_selection = df_selection.astype("object")
        self.dfs["selection"] = df_selection.set_index("selection_id")
        
        # Table of all selections defined in each batch
        key_order = ["batch", "method", "dataset", "selection_id", "data_id", "file_names"]
        df = pd.DataFrame(self._get_combined_configurations(dataset_param_combs, selection_param_combs))
        df["file_names"] = df.apply(
            lambda x: f"selection/{x['method']}_{x['selection_id']}_{x['dataset']}_{x['data_id']}.csv", axis=1
        )
        for key in key_order[::-1]:
            df.insert(0, key, df.pop(key))
        df = df.astype("object")
        self.dfs["selection_overview"] = df
        
        # Save tables
        self.dfs["data"].to_csv(self.data_params_file)
        self.dfs["selection"].to_csv(self.selection_params_file)
        self.dfs["selection_overview"].to_csv(self.selection_overview_file)
        
        # Get all pipeline output file names
        self.file_names = self.dfs["selection_overview"]["file_names"].unique().tolist()
        #self.file_names = [str(Path(self.RESULTS_DIR, file_name)) for file_name in file_names]
        
   
    def get_selection_params(self, selection_id: int) -> dict:
        """Get the selection parameters for a given selection id
        
        Arguments
        ---------
        selection_id: int
            The selection id
        
        Returns
        -------
        selection_params: dict
            The selection parameters
        """
        selection_params = self.dfs["selection"].loc[selection_id].to_dict()
        return selection_params
    
    def get_data_params(self, data_id: int) -> dict:
        """Get the dataset parameters for a given dataset id
        
        Arguments
        ---------
        data_id: int
            The dataset id
        
        Returns
        -------
        data_params: dict
            The dataset parameters
        """
        data_params = self.dfs["data"].loc[data_id].to_dict()
        return data_params
        
    def _get_hyperparams_that_occur_in_config(self) -> dict:
        """
        
        Note: this function is not used anymore since we just use all default hyperparameters
        
        Returns
        -------
        dictionary like DEFAULT_PARAMETERS but only with hyperparameters that occur in the config
        }
        """
        hyperparams = {"dataset":[], "selection":[]}
        for batch, batch_dict in self.config["selections"].items():
            if "dataset_param" in batch_dict.keys():
                for param in batch_dict["dataset_param"].keys():
                    if param not in DEFAULT_PARAMETERS["dataset"].keys():
                        raise ValueError(f"Parameter {param} ({batch}, dataset) not in DEFAULT_PARAMETERS")
                    if param not in hyperparams["dataset"]:
                        hyperparams["dataset"].append(param)
            if "selection_param" in batch_dict.keys():
                for param in batch_dict["selection_param"].keys():
                    if param not in DEFAULT_PARAMETERS["selection"].keys():
                        raise ValueError(f"Parameter {param} ({batch}, selection) not in DEFAULT_PARAMETERS")
                    if param not in hyperparams["selection"]:
                        hyperparams["selection"].append(param)
        
        default_hparams = DEFAULT_PARAMETERS.copy()
        for key in default_hparams.keys():
            default_hparams[key] = {k:v for k,v in default_hparams[key].items() if k in hyperparams[key]}
                                
        return default_hparams
    
    
    def _get_combinations_of_params(
            self, 
            batch_dict: dict, 
            default_hparams: dict,
            main_key : str = "methods",
            main_key_param_name : str = "method",
            param_key : str = "selection_param",
            default_param_key : str = "selection",
        ) -> List[dict]:
        """Convert config into list of dictionaries with all combinations of parameters
        
        example how the config looks like:
        batch2:
            ...
            selection_param:
                n: [50,100,150]
                penalty : [None, "highly_expressed_penalty"]
            methods:
                spapros
                nsforest
                
        Arguments
        ---------
        --> batch_dict = {
                "selection_param" : {
                    "n": [50,100,150],
                    "penalty" : [None, "highly_expressed_penalty"]
                },
                "methods": ["spapros", "nsforest"]
            }
        
        Returns
        -------
        --> param_dicts = [
                {'method': 'spapros', 'n': 50, 'penalty': None},
                {'method': 'spapros', 'n': 50, 'penalty': 'highly_expressed_penalty'},
                ...
            ]

        """
        
        mkp_name = main_key_param_name
        if isinstance(batch_dict[main_key], list):
            names = batch_dict[main_key]
        elif " " in batch_dict[main_key]:
            names = batch_dict[main_key].split(" ")
        else:
            names = [batch_dict[main_key]]
        
        if param_key not in batch_dict.keys():
            param_dicts = [{mkp_name:n} for n in names]
        else:
            # Get all possible combinations of dataset parameters as list of dictionaries
            val_combs = list(itertools.product(*batch_dict[param_key].values()))
            param_dicts = [{key:val for key, val in zip(batch_dict[param_key].keys(), vals)} for vals in val_combs]
            
            # Repeat combinations for each of the datasets and add the dataset name
            param_dicts = [
                {mkp_name:n, **p_dict} for n, p_dict in list(itertools.product(names, param_dicts))
            ]
            
        # Add default hyperparameters that are not specified in the config
        for param_dict in param_dicts:
            for key, val in default_hparams[default_param_key].items():
                if key not in param_dict.keys():
                    param_dict[key] = val
                    
        # Order each dictionary by the order of [mkp_name] + default_hparams[default_param_key].keys()
        param_dicts = [
            {key: param_dict[key] for key in [mkp_name] + list(default_hparams[default_param_key].keys())} for param_dict in param_dicts
        ]
        
        return param_dicts    
    
    
    def _convert_pseudo_params_to_param(
            self, 
            config: dict,

        ) -> List[dict]:
        """Convert pseudo parameters to actual parameters
        
        Arguments
        ---------
        config: dict
            The config dictionary as read from the yaml file
            
        Returns
        -------
        dict
            The config dictionary with the pseudo parameters converted to actual parameters
        
        """
        # Get pseudo parameters
        pseudo_params = {
            "dataset" : [v for v in PSEUDO_PARAM_TO_PARAM["dataset"].keys()],
            "selection" : [v for v in PSEUDO_PARAM_TO_PARAM["selection"].keys()],
        }
        
        # Get all parameters that are set via pseudo parameters
        target_params = {
            "dataset" : [v for v in PSEUDO_PARAM_TO_PARAM["dataset"].values()],
            "selection" : [v for v in PSEUDO_PARAM_TO_PARAM["selection"].values()],
        }
        
        new_config_s = {}
        old_config_s = config["selections"]
        # Iterate over each batch of the configuration
        for batch in old_config_s:
            new_config_s[batch] = {}
            
            # Copy datasets and methods
            new_config_s[batch]["datasets"] = old_config_s[batch]["datasets"]
            new_config_s[batch]["methods"] = old_config_s[batch]["methods"]
            
            # Update or copy dataset parameters
            if "dataset_param" in old_config_s[batch]:
                new_config_s[batch]["dataset_param"] = {}
                
                for param in old_config_s[batch]["dataset_param"]:
                    # Raise an error if a parameter is given that should be set via its pseudo parameter
                    if param in target_params["dataset"]:
                        raise ValueError(
                            f"Parameter {param} ({batch}, dataset) can not be set directly, \
                            it is set via the pseudo parameter {PARAM_TO_PSEUDO_PARAM['dataset'][param]}"
                        )
                    # Convert pseudo parameters
                    elif param in pseudo_params["dataset"]:
                        new_config_s[batch]["dataset_param"].update({
                            PSEUDO_PARAM_TO_PARAM["dataset"][param] : self._dataset_pseudo_parameter_conversion(
                                param, old_config_s[batch]["dataset_param"][param]
                            )
                        })
                    # Copy parameters that are not converted
                    else:
                        new_config_s[batch]["dataset_param"].update({
                            param : old_config_s[batch]["dataset_param"][param]
                        })
            
            # Update or copy selection parameters
            if "selection_param" in config["selections"][batch]:
                new_config_s[batch]["selection_param"] = {}
                
                for param in old_config_s[batch]["selection_param"]:
                    # Raise an error if a parameter is given that should be set via its pseudo parameter
                    if param in target_params["selection"]:
                        raise ValueError(
                            f"Parameter {param} ({batch}, selection) can not be set directly, \
                            it is set via the pseudo parameter {PARAM_TO_PSEUDO_PARAM['selection'][param]}"
                        )
                    # Convert pseudo parameters
                    elif param in pseudo_params["selection"]:
                        new_config_s[batch]["selection_param"].update({
                            PSEUDO_PARAM_TO_PARAM["selection"][param] : self._selection_pseudo_parameter_conversion(
                                param, old_config_s[batch]["selection_param"][param]
                            )
                        })
                    # Copy parameters that are not converted
                    else:
                        new_config_s[batch]["selection_param"].update({
                            param : old_config_s[batch]["selection_param"][param]
                        })
                        
        # Update the config dictionary
        config["selections"] = new_config_s
        
        return config
                        
        
    def _dataset_pseudo_parameter_conversion(self, param: str, values: list) -> list:
        """Convert pseudo parameters to actual parameters
        
        Arguments
        ---------
        param: str
            The pseudo parameter to convert
        values: list
            The values of the parameter to convert
            
        Returns
        -------
        list
            The converted values
        
        """
        if param == "cells_per_ct_n_seeds":
            return [seed for seed in range(values[0])]

    def _selection_pseudo_parameter_conversion(self, param: str, values: list) -> list:
        """Convert pseudo parameters to actual parameters
        
        Arguments
        ---------
        param: str
            The parameter to convert
        values: list
            The values of the parameter to convert
            
        Returns
        -------
        list
            The converted values
        
        """
        # So far no pseudo parameters for selection parameters
        pass
        
        
    def _get_ids_and_configs(
            self, 
            param_combs: List[dict], 
            id_str: str = "id", 
            ids_to_config_old: Optional[dict] = None
        ) -> Tuple[dict, dict]:
        """
        
        Returns
        -------
        ids_to_config: dict
            Dictionary with ids as keys and the respective configuration as values
            e.g.: {
                0 : {'dataset': 'sc_mouse_brain', 'n_cts': None, 'cells_per_ct': None}
                1 : {'dataset': 'sc_mouse_brain', 'n_cts': None, 'cells_per_ct': 50}
            }
        param_combs: dict
            Dictionary with batch names as keys and a list of dictionaries with all combinations of parameters as 
            values. Now also including the configuration id.
            e.g.: {
                "batch1" : [
                    {'dataset': 'sc_mouse_brain', 'n_cts': None, 'cells_per_ct': None, 'id': 0},
                    {'dataset': 'sc_mouse_brain', 'n_cts': None, 'cells_per_ct': 50, 'id': 1},
                ],
                "batch2" : [
                    {'dataset': 'sc_mouse_brain', 'n_cts': None, 'cells_per_ct': None, 'id': 0}
                ]
        
        """
        p_combs = param_combs.copy()
        
        ids_to_config = ids_to_config_old if ids_to_config_old is not None else {}
        idx = 0 if ids_to_config_old is None else max(ids_to_config_old.keys()) + 1
        for batch, configs in p_combs.items():
            for i, config in enumerate(configs):
                if config not in ids_to_config.values():
                    ids_to_config[idx] = config.copy()
                    p_combs[batch][i][id_str] = idx
                    idx += 1
                else:
                    p_combs[batch][i][id_str] = list(ids_to_config.keys())[
                        list(ids_to_config.values()).index(config)
                    ]        
        
        return ids_to_config, p_combs
      
      
    def _load_config_table(self, config_file: str, param_group: str = "selection") -> dict:
        """Load the config table and convert it to a dictionary
        
        Arguments
        ---------
        config_file: str
            The path to the config file
        param_group: str
            Either "dataset" or "selection"
        
        Returns
        -------
        config_dict: dict
            Dictionary with ids as keys and the respective configuration as values
        
        """
        df = pd.read_csv(config_file, index_col=0)
        df = df.astype("object")
        for col in df.columns:
            df.loc[df[col].isnull(), col] = "None"
            df.loc[df[col] != "None", col] = df.loc[df[col] != "None", col].astype(DEFAULT_PARAMETERS_TYPES[param_group][col]).tolist()
            
        config_dict = df.to_dict(orient='index')
        
        return config_dict
        
    def _get_combined_configurations(self, dataset_param_combs: dict, selection_param_combs: dict) -> List[dict]:
        """Get all combinations of dataset and selection configurations within each batch
        """
        
        combined_configs = {}
        for batch in dataset_param_combs.keys():
            dataset_configs = dataset_param_combs[batch]
            selection_configs = selection_param_combs[batch]
            
            combs = list(itertools.product(dataset_configs, selection_configs))
            combined_configs[batch] = [
                {**dataset_config, **selection_config, **{"batch":batch}} for dataset_config, selection_config in combs
            ]
        
        combined_configs = [
            config for batch in combined_configs.keys() for config in combined_configs[batch]
        ]
        
        return combined_configs
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    #def _get_combinations_of_dataset_params(self, batch_dict: dict, default_hparams: dict) -> List[dict]:
    #    """Convert dataset config into list of dictionaries with all combinations of dataset parameters
    #    
    #    example how the config looks like:
    #    batch2:
    #        ...
    #        datasets: [sc_mouse_brain, sc_human_brain]
    #        dataset_param:
    #            cells_per_ct: [50,100,200,500,1000,None]
    #            gene_key: ["hvg_probe_constraint"]
    #    
    #    Arguments
    #    ---------   
    #    --> batch_dict = {  
    #            "datasets"      : ["sc_mouse_brain", "sc_human_brain"],
    #            "dataset_param" : {
    #                "cells_per_ct": [50,100,200,500,1000,None],
    #                "gene_key": ["hvg_probe_constraint"]
    #            }
    #        }
    #    
    #    Returns
    #    -------
    #    --> param_dicts = [
    #            {'datasets': 'sc_mouse_brain', 'cells_per_ct': 50, 'gene_key': 'hvg_probe_constraint'},
    #            {'datasets': 'sc_mouse_brain', 'cells_per_ct': 100,'gene_key': 'hvg_probe_constraint'},
    #            ...
    #        ]
    #    """
    #    
    #    if "dataset_param" not in batch_dict.keys():
    #        param_dicts = [{"dataset":dataset} for dataset in batch_dict['datasets']]
    #    else:
    #        # Get all possible combinations of dataset parameters as list of dictionaries
    #        val_combs = list(itertools.product(*batch_dict['dataset_param'].values()))
    #        param_dicts = [{key:val for key, val in zip(batch_dict['dataset_param'].keys(), vals)} for vals in val_combs]
    #        
    #        # Repeat combinations for each of the datasets and add the dataset name
    #        param_dicts = [
    #            {"dataset":dataset, **p_dict} for dataset, p_dict in list(itertools.product(batch_dict['datasets'], param_dicts))
    #        ]
    #        
    #    # Add default hyperparameters that are not specified in the config
    #    for param_dict in param_dicts:
    #        for key, val in default_hparams["dataset"].items():
    #            if key not in param_dict.keys():
    #                param_dict[key] = val
    #                
    #    # Order each dictionary by the order of ["dataset"] + default_hparams["dataset"].keys()
    #    param_dicts = [
    #        {key: param_dict[key] for key in ["dataset"] + list(default_hparams["dataset"].keys())} for param_dict in param_dicts
    #    ]
    #    
    #    
    #    return param_dicts
    #
    #def _get_combinations_of_selection_params(self, batch_dict: dict, default_hparams: dict) -> List[dict]:
    #    """Convert selection config into list of dictionaries with all combinations of selection parameters
    #    
    #    example how the config looks like:
    #    batch2:
    #        ...
    #        selection_param:
    #            n: [50,100,150]
    #            penalty : [None, "highly_expressed_penalty"]
    #        methods:
    #            spapros
    #            nsforest
    #            
    #    Arguments
    #    ---------
    #    --> batch_dict = {
    #            "selection_param" : {
    #                "n": [50,100,150],
    #                "penalty" : [None, "highly_expressed_penalty"]
    #            },
    #            "methods": ["spapros", "nsforest"]
    #        }
    #    
    #    Returns
    #    -------
    #    --> param_dicts = [
    #            {'method': 'spapros', 'n': 50, 'penalty': None},
    #            {'method': 'spapros', 'n': 50, 'penalty': 'highly_expressed_penalty'},
    #            ...
    #        ]
    #        
    #    """
    #    
    #    if "selection_param" not in batch_dict.keys():
    #        param_dicts = [{"method":method} for method in batch_dict['methods']]
    #    else:
    #        # Get all possible combinations of selection parameters as list of dictionaries
    #        val_combs = list(itertools.product(*batch_dict['selection_param'].values()))
    #        param_dicts = [{key:val for key, val in zip(batch_dict['selection_param'].keys(), vals)} for vals in val_combs]
    #        
    #        # Repeat combinations for each of the methods and add the method name
    #        param_dicts = [
    #            {"method":method, **p_dict} for method, p_dict in list(itertools.product(batch_dict['methods'], param_dicts))
    #        ]
    #        
    #    # Add default hyperparameters that are not specified in the config
    #    for param_dict in param_dicts:
    #        for key, val in default_hparams["selection"].items():
    #            if key not in param_dict.keys():
    #                param_dict[key] = val
    #    
    #    # Order each dictionary by the order of ["method"] + default_hparams["selection"].keys()
    #    param_dicts = [
    #        {key: param_dict[key] for key in ["method"] + list(default_hparams["selection"].keys())} for param_dict in param_dicts
    #    ]
    #    
    #    return param_dicts
