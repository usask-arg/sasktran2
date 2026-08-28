
(_api_database)=
# Databases

## Database Storage

By default, SASKTRAN2 stores downloaded and generated database files in the
platform-specific user data directory. Set the `SASKTRAN2_DATABASE_ROOT`
environment variable to use a different directory:

```bash
export SASKTRAN2_DATABASE_ROOT=/opt/sasktran2/database
```

The environment variable takes precedence over the `database_root` value in the
user configuration file. An explicit `db_root` passed to an individual database
object takes precedence over both global settings.

## Web Databases
```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.database.StandardDatabase
```

## Internal Databases
```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.database.MieDatabase
```

## Solar Databases
```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.solar.SolarModel
```
