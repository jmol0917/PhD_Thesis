# GEMNAST_pub
Computational high throughput platform for the assessment of nutritional and metabolic attributes of microorganisms/cells as represented in genome scale models 

The scripts in this repository are required to replicate the results from the Thesis titled "" and include the corresponding reactions and dictionaries.

The scripts are ready to be ran once the correct directory inputs and outputs are specified. Further, download of the corresponding GSM models in SBML format 
from www.vmh.life is required and the directory where these files are located should be provided as an input.

To download AGORA SBML files:
1. Go to https://www.vmh.life/#home
2. Click on 'Gut Microbiota'
3. Click on 'Download' at the bar near the TOP of the webpage
4. Click on 'Microbe collections'
5. Click on the 'AGORA 1.03 without mucin' version
6. All the models within this version will be downloaded as part of a .zip file

GEMNAST is divided in two branches. The Non-Replaceable Nutrients Assessment and the Replaceable Nutrients Assessment.

The Non-Replaceable Nutrients Assessment is composed of three scripts:
The essential_nutrients_assessment.py script produces individual files for every examined GSM. The essential_nutrients_compilation.py script can be used 
to generate a single boolean table with information from all the resulting files.

The optional_nutrients_assessment.py and the optional_nutrients_export.py files already produce a single boolean table with information from all the assessed models and
therefore additional scripts are not required.

The Replaceable Nutrients Assessment is composed of three scripts:
The replaceable_nutrients_degradation_assessment.py script produces...

The repalceable_nutrient_degradation_and_intermediate_metabolite_synthesis.py script produces...

The repalceable_nutrient_degradation_and_intermediate_metabolite_synthesis.py script produces...

For questions or comments please email juan.molina@sydney.edu.au
