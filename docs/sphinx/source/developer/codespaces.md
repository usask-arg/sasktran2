(_dev_codespaces)=
# Using GitHub Codespaces
An easy way to contribute small changes to SASKTRAN2 is through GitHub codespaces.
This is suitable for fixing small typos or easy to spot bugs, adding new tests, or improving the documentation.
For more in-depth development work it is recommended to set up your own local development environment.

## Creating a Codespace
When working with a codespace it is recommended to do development on your own local fork. If you haven't already forked
the repository, you can go to https://github.com/usask-arg/sasktran2 and press the fork button on the top right.  This will create a copy of the repository in your local Github namespace.

Next, go to your forked repository.  The URL will look something like https://github.com/USERNAME/sasktran2, on the top right click the green
"Code" button and press the codespace tab and then create a new codespace.

Once inside the codespace, go to the terminal at the bottom and run

```
micromamba activate sasktran2-dev-env
pip install -e .
```

which will compile the code.  Note that it is necessary to re-run `pip install -e .` everytime you change any of the c++ code inside the model.

## Running the Python Tests
First, setup vscode to use the correct python environment.  Run `ctrl+shift+p`, and select `Python: select interpreter` and choose `sasktran2-dev-env`.

Next, on the left, click on the `testing` button and press configure Python tests.  Choose `pytest` as the testing framework and choose the `tests` folder
when prompted.  You can then press the run button next to `sasktran2` to run all of the tests.

## Building and Viewing the Documentation
First, go to the extensions tab on the left, and install the "Live Preview" extension.

Then, in the terminal, run

```
cd docs/sphinx
make html
```

Then, on the explorer on the left, go to `docs/sphinx/build` and right click on `index.html` and press "Show Preview". A button
should appear prompting you to open a web browser.  On this page you may have to navigate to `docs/sphinx/build` again.

## Making Changes
To do development inside the codespace, press the "Main" button on the bottom left, and choose "Create new branch".  Here you can make
a new branch.  Then once you are done making changes, press the "Source Control" button on the left and commit your changes.  Make sure to
publish the changes to your fork after committing.  Then you can make a pull request as normal.
