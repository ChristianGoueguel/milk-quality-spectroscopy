# Milk Diagnostics Dataset 

This is the dataset of sensor spectra coupled with lab results for milk diagnostics. The dataset is organized by sensor and includes dark spectra
-  samples
-  and lab results. Each sample is associated with a specific milk session and contains both raw spectra and lab results. 

## Organization

The dataset is organized as follows:
 * dataset/sensor_{stall_id}
 * dataset/sensor_{stall_id}/lab_results.csv **tube_no is the primary key for the lab results and relates to the tube_no in the spectra.csv files**
 * dataset/sensor_{stall_id}/sensor.csv
 * dataset/sensor_{stall_id}/dark_spectra/dark_spectra.parquet. (zstd compression)
 * dataset/sensor_{stall_id}/spectra/tube_no_{tube_no}.parquet. (zstd compression)



## GOTCHAS about this data

* Temperature on each sensor is not standardized across sensors. Temperature is not in a standard unit, is a raw ADC value.
* LED Lamp current measurement is not standardized across sensors, any value shift here must be considered only in relation to this sensor.  It is a raw ADC value.
* Each sensor has its own set of wavelenghts for all its measured spectra.

* ~~!!! Some Lab Results do not yet have corresponding spectra. This is a known issue and will be addressed in future updates to this dataset.~~
* ~~Correlary to the above ^^ Some milkings have more than one animal in the spectra data. This is a known issue and will be addressed in future updates to this dataset.~~
* ~~How to spot the problem above:~~
  - ~~If the tube_no in the lab_results.csv does not have a corresponding file in the spectra directory, it is a problem.~~
  - ~~If the spectra file has more than about 13 minutes of data it likely contains two animals worth, and needs to be split differently and reorganized.~~
  - ~~I believe there are about 12 or 13 instances in this dataset~~
* I have implemented a 'reselection' that uses the timestamp from the milking parlor as to when the cow was removed from the milker, MOSTLY this works well.
* ~~I've noticed a few situatinos where some cows have a huge milk output but we don't see it, I'm not sure what to make of this, if we just never saw the milk because of some nyquist sampling issue or if it's a real issue in the data that nees to be addressed.~~ 
* ^^ I found a that the clock sources were off by a bit of time, and we're also probably missing times when the tube is full of milk.  Such is the data. 
* ^^^ I found that the Endtime timestamp is actually a Starttime, and that sorted out a LOT of the weird data that I was seeing.

## CHANGELOG

1. **Initial Release**
   - Initial release of the dataset.
1. **Full Dataset**
   - Full dataset with all the data.
   - The spectra are now in parquet format
   - nm_1234* columns have been replaced with a 'spec_array' array column to save having to pivot data, the corresponding wavelengths are in the 'wavelengths' column in sensor.csv
1. ** Updates **
   - Removed the "most empty" spectra from the dataset, there's no milk in the tube, there's nothing to find in there.
   - Found some missing data and added it back.
   - Adjusted timestamp selection to recognize that "endtime" timestamp is actually a starttime.
   - Stall 4 nearly never clears the milk tube, so I am curious to see the differences in that vs the other sensors


## Files

  * dataset/sensor_{stall_id}/lab_results.csv
    - Contains the lab results for each sample.
    - tube_no is the primary key for the lab results and relates to the tube_no in the spectra.csv files.
    - Columns are:
      - tube_no # tube number written on the tube for the lab
      - barnname # Cow Ear Tage
      - rfid # RFID Tag for cow
      - lact # number of times cow has been successfully bred to milk
      - grp # group
      - dim_days # days since the cow started producing milk
      - date # date of the sample
      - day_of_study # day of the study starting at 1
      - daily_milking # Number of times this cow has been milked today 1, 2, and sometimes 3
      - study_milking # Number of times this cow has been milked during the study
      - milking_end # End time of the milking session
      - avgmilkdur_mins # Average milking duration in minutes, according to the milk weight system
      - milking_stall # which stall the cow was in during the milking session
      - milk_wt_lbs # weight of the milk in pounds
      - tube_number # tube number for the lab sample, same as tube_no
      - date_1 # date of the sample, same as date
      - fat_percent # fat percentage in the milk as measured by the lab
      - protein_percent # protein percentage in the milk as measured by the lab
      - scc_thous_per_ml # somatic cell count in thousands per milliliter as measured by the lab
      - lactose_percent # lactose percentage in the milk as measured by the lab
      - other_solids_percent # other solids percentage in the milk as measured by the lab
      - mun_mg_per_dl # total solids in milligrams per deciliter as measured by the lab
      - solids_not_fat_percent # total solids that are not fat as a percentage in the milk as measured by the lab
      - total_solids_percent # total solids as a percentage in the milk as measured by the lab
      - endtime # datetime derrived from date and milking_end
  * dataset/sensor_{stall_id}/sensor.csv
    - Contains information about the sensor including the specific wavelenghts and other details
    - columns:
      - sensor_uid # the unique identifier for the sensor
      - pcba # the printed circuit board assembly number for the sensor
      - box_number # the outer assembly box number for the sensor
      - adc_gain # the analog-to-digital converter gain for the VIS-NIR imager
      - adc_offset # the analog-to-digital converter offset for the VIS-NIR imager, centered around 512, anything less than 512 is a negative offset, anything greater than 512 is a positive offset, 511 being the greatest negative offset, 1024 being the greatest positive offset
      - stall_number # the stall number where the sensor is located, same as lab_results.csv.milking_stall
      - stall_title # Human readable stall title
      - imager_sn # the serial number of the VIS-NIR imager 
      - sensor_wavelength_coefficients # the wavelength coefficients for the sensor, used to derrive the actual measured wavelengths this sensor sees in nanometers
      - wavelengths_nm # the actual measured wavelengths this sensor sees in nanometers, as derrived from sensor_wavelength_coefficients
  * dataset/sensor_{stall_id}/spectra/tube_no_{tube_no}.csv
    - Contains all the spectra that is coupled with this lab sample
    - The `spec_type` column indicates the type of spectrum, dark, sample or empty.  This is a a rough observation based on the data read from the sensor and is not an authorative source.
    - Columns:
      - tube_no # tube_number, same as lab_results.csv
      - old_spec_id # spectrum_id for this milking_session, because of the data errors described above, this number will not be unique and will not be stable for a give spectrum
      - datetime # the datetime that the spectrum was recorded
      - spec_type # the type of spectrum, dark, sample or empty, NOT AUTHORATIVE
      - session_from # Either 'reselected_from_timestamps' or 'sensor_detected' generally the 'sensor_detected' are more consistent with each other.
      - temp # temperature, raw data from the sensor, not comparable between sensors today
      - integration_time # the integration time used for this measurement, in number of clocks of some defined clock in the sensor. 250 is somewhat less than a tenth of a second
      - led_current_setting_mA # the target milliamp current for the LED Lamp source
      - led_current_measured_raw # raw analog to digital reading of the LED lamp current, each sensor is mildly different and so this reading is NOT compariable across sensors
      - spec_array  # raw 16 bit integer value array from eash wavelength from the sensor.

  * dataset/sensor_{stall_id}/dark_spectra/dark_spectra.csv
    - Contains all the dark spectra as taken by this sensor.  The schema is the same as the spectra/tube_no_{tube_no}.csv files, but there is no `tube_no` column

 
