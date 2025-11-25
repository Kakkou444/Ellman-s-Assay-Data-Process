**Plate Reader Data Processing Pipeline**  
This repository contains two Python scripts for processing Gen5 plate reader data and generating dose–response curves for inhibitor assays.
1. bulkprocess.py:
  Process raw data (.xlsx or .csv file) obtained from Gen5 plate reader. Create a folder "bulkdata" in your working folder, throw all your data in the "bulkdata", the code will process the data to give absorbance change rate and calculate %activity and %inhibitor, and save as .csv file with the corresponding file name.
2. dose_resp_curve.py:
  Process the .csv data produced from bulkprocess.py, generate a summary with mean %activity for 3 replicates and the SD, and produce a dose response curve using 4PL fitting. The program requires you to read the file manually (change the file name manually at the beginning of the code), but it would save the file with the corresponding name automatically.

P.S. The program can only work with column 1 on the 96 plate set as no inhibitor control and column 2 as no enzyme control. Please delete any extra information in your .xlsx or .csv file before feeding it into bulkprocess.py. Example plate layout and raw data format are shown in the graph.
<img width="1325" height="578" alt="image" src="https://github.com/user-attachments/assets/989be82b-f407-4f43-ace2-8234051db2f1" />
<img width="2493" height="956" alt="image" src="https://github.com/user-attachments/assets/452dab8f-e437-4904-b062-5a2a416eb968" />
