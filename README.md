# ModuleCalibration
This repository contains the code used to create the 2D maps and 1D histograms starting from the root file which are the direct output of Ph2ACF.

All the scripts takes all the input parameters from the `default_config.yaml`. There is also the possibility to specify your custom parameters inside another `custom_config.yaml` that you should parse in the command line. For example, if you want to run the script `single_modules.py` with a custom configuration you can run the following command.

```
python pixelalive.py custom_config.yaml
```

In this way the custom and the default yaml files will be merged together, and if the same parameters are specified in both configs the `custom_config.yaml` will overwrite the default. For this reason it is recommended to keep the `default_config.yaml` unchanged and to define your own `custom_config.yaml`.

Inside the `custom_config.yaml` you will have to specify
* your input folder and root files:
  ```
  input:
    input_folder: /home/giorgia/module_calibration/
    input_file: scurve.root
  ```
  note that the root file given has to be the output of Ph2ACF for that specific scan.
  
* your output folder:
  ```
  output:
    output_folder: /home/giorgia/module_calibration/plots/
  ```
* the chip IDs of your tested module, included if it is dual or quad, and the hybrid to which it is connected:
  ```
  ph2acf:
    chip_id: [0,1]
    hybrid: 1
    is_dual: False
  ```    
