# Figures from Petersen, Roule, Fouvry, Pichon, Tep (2024)

This is all of the python code necessary to generate the figures, including the data files[^1]. The scripts themselves are located in `scripts/`, while pre-generated versions of the figures are in `figures/`. Figure numbers follow the text. Inside of `scipts/`, another directory `data/` holds the files needed to generate the figures, divided up by the figure numbers.

Some scripts will require the external `exptool` library in order to inferface with the EXP simulation outputs, available here: https://github.com/michael-petersen/exptool.

For data files too large to be included, we provide the scripts that generated the data files in the `scripts/data/datagenerators/` directory. Other figure generating scripts may be found in the julia libraries themselves.


[^1]: Except for some data files that are too big for GitHub; the user will need to generate themselves. See below.
