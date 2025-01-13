from __future__ import annotations

import logging
import urllib.request
import zipfile
from pathlib import Path

import appdirs
import yaml
from packaging import version

APPDIRS = appdirs.AppDirs(appname="sasktran2", appauthor="usask-arg")


def user_config_file_location() -> Path:
    """
    File location of the user config file
    """
    user_config_dir = APPDIRS.user_config_dir

    return Path(user_config_dir).joinpath("config.yml")


def load_user_config():
    """
    Opens the user config file as a dictionary
    """
    user_config_file = user_config_file_location()

    try:
        with Path.open(user_config_file) as f:
            if version.parse(yaml.__version__) > version.parse("5"):
                config = yaml.load(f, Loader=yaml.FullLoader)
            else:
                config = yaml.load(f)

        if config is not None:
            return config
        return {}
    except FileNotFoundError:
        return {}


def save_user_config(user_config: dict):
    """
    Saves the user config from a dictionry
    """
    user_config_file = user_config_file_location()

    Path(user_config_file).parent.mkdir(exist_ok=True, parents=True)

    with Path.open(user_config_file, "w") as f:
        yaml.dump(user_config, f, default_flow_style=False)


def database_root() -> Path:
    dir = load_user_config().get("database_root", None)

    if dir is None:
        return Path(APPDIRS.user_data_dir).joinpath("database")
    return Path(dir)


def are_extended_db_downloaded() -> bool:
    """
    Checks if the extended databases are downloaded
    """
    if database_root() is None:
        return False
    return database_root().joinpath("cross_sections/h2o/hitran_01nm_res.nc").exists()


def download_extended_databases(version: str = "latest"):
    db_root = database_root()

    if db_root is None:
        msg = "Database root not set"
        raise OSError(msg)

    output_file_temp = db_root.joinpath(f"v_{version}_extended.zip")

    if output_file_temp.exists():
        output_file_temp.unlink()

    try:
        urllib.request.urlretrieve(
            f"https://arg.usask.ca/sasktranfiles/sasktran2_db/v_{version}_extended.zip",
            filename=output_file_temp.as_posix(),
        )
        with zipfile.ZipFile(output_file_temp.as_posix(), "r") as zip_ref:
            zip_ref.extractall(db_root)
    except Exception as e:
        logging.exception(e)

    if output_file_temp.exists():
        output_file_temp.unlink()
