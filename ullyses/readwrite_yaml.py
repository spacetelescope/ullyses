import yaml
import glob
import os

def read_config(yamlfile):
    """
    Read a yaml configuration file for STIS re-calibration parameters.
    Args:
        target (str): Name of target, which will determine yaml file to use.
        yamlfile (str): Name of yaml file to read in.
    Returns:
        data (dict): Dictionary with all custom calibration pars.
    """

    with open(yamlfile) as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
        
        return data

def write_config(yaml_d, yamlfile):
    """
    Write final yaml configuration file for STIS re-calibration parameters.
    Args:
        yaml_d (dict): Dictionary to be written as a YAML file.
        yamlfile (str): Name of yaml file to write out.
    Returns:
        None
    """

    with open(yamlfile, "w") as f:
        data = yaml.dump(yaml_d, f)
