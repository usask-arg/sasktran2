import logging
import urllib.request
import zipfile
from pathlib import Path

import appdirs
import yaml
from packaging import version

APPDIRS = appdirs.AppDirs(appname='sasktran2', appauthor='usask-arg')

def user_config_file_location() -> Path:
    """
    File location of the user config file
    """
    user_config_dir = APPDIRS.user_config_dir

    return Path(user_config_dir).joinpath('config.yml')



def load_user_config():
    """
    Opens the user config file as a dictionary
    """
    user_config_file = user_config_file_location()

    try:
        with open(user_config_file) as f:
            if version.parse(yaml.__version__) > version.parse('5'):
                config = yaml.load(f, Loader=yaml.FullLoader)
            else:
                config = yaml.load(f)

        if config is not None:
            return config
        else:
            return {}
    except FileNotFoundError:
        return {}


def save_user_config(user_config: dict):
    """
    Saves the user config from a dictionry
    """
    user_config_file = user_config_file_location()

    Path(user_config_file).parent.mkdir(exist_ok=True, parents=True)

    with open(user_config_file, 'w') as f:
        yaml.dump(user_config, f, default_flow_style=False)


def database_root() -> Path:
    dir = load_user_config().get('database_root', None)

    if dir is None:
        return None
    else:
        return Path(dir)


def download_standard_databases(version : str = 'latest'):
    user_config = load_user_config()

    if 'database_root' in user_config:
        data_directory = Path(user_config['database_root'])
    else:
        data_directory = Path(APPDIRS.user_data_dir).joinpath('database')

    if not data_directory.exists():
        data_directory.mkdir(parents=True)

    output_file_temp = data_directory.joinpath(f'v_{version}.zip')

    if output_file_temp.exists():
        output_file_temp.unlink()

    try:
        urllib.request.urlretrieve('https://arg.usask.ca/sasktranfiles/sasktran2_db/v_{}.zip'.format(version),
                                    filename=output_file_temp.as_posix())
        with zipfile.ZipFile(output_file_temp.as_posix(), 'r') as zip_ref:
            zip_ref.extractall(data_directory)
    except Exception as e:
        logging.exception(e)

    if output_file_temp.exists():
        output_file_temp.unlink()

    user_config['database_root'] = data_directory.as_posix()
    save_user_config(user_config)

if __name__ == "__main__":
    download_standard_databases()