# IREDFisher: A web server for fast prediction of imine reductases hits with given substrates.

Chiral amines are common building blocks for producing pharmaceuticals and commodity chemicals. Imine reductases (IRED) show attractive properties in the broad prochiral substrates scope and has become one of the hottest biocatalyst in chiral amine production. However, finding hits by large-scale screening is inefficient and cost much in time and money. Here we propose our computational workflow IREDFisher to optimize the IRED screening panel with the analysis on sequence and structure. Tests based on published data showed the workflow’s capability to enrich the most active enzymes in top 20 ranks and to improve the screening efficiency for finding hits with high conversion by up to 3 folds, which is further confirmed by our lab validations on four different substrates.  A web interface was built and provided to the public free of charge. With proved effectiveness and easy inputs, IREDFisher is expected to benefit the biocatalysis community.
</br>
</br>
![image](https://user-images.githubusercontent.com/61114239/116463260-659d4880-a862-11eb-88f8-3ce8620f98aa.png) 
##### Figure 1. Flow chart of the four-step IREDFisher workflow



![image](https://user-images.githubusercontent.com/61114239/116463288-6fbf4700-a862-11eb-8644-1660dcf8c018.png)
##### Figure 2. Web interface of IREDFisher. (a) Job submission interface (b) Input file examples. Top panel: fasta-formatted sequence file. Bottom panel: generation of substrate structure and N atom index by Marvin JS (https://marvinjs-demo.chemaxon.com/latest/) (c) Output files. The prediction sequence name file (1a_imine_sort_rescore.csv in the example) is shown.
 </br>


As there are many softwares needed for running IREDFisher, we recommend to use our webserver at https://enzymeevolver.com/IREDFisher 
or run our built docker container on your local machine.

## Input file preparation
1. sequence file: fasta format file </br>
  ![image](https://user-images.githubusercontent.com/61114239/116700404-a6aa6f80-a9be-11eb-994d-d953a2e323bd.png)
2. Ligand preparation. </br>
   The imine form should be used as the input structure. <br>
   i. go to Marvin JS online molecule drawer https://marvinjs-demo.chemaxon.com/latest/demo.html </br>
   ii. draw the imine substrate.</br>
   iii. add hydrogens and display the atom index </br>
   iv. save as pdb formatted file.</br>
   ![image](https://user-images.githubusercontent.com/61114239/116701229-934bd400-a9bf-11eb-86cc-919a3f12391b.png)
   


## Run IREDFisher Docker container
Note:
Please make sure that you have got license from MODELLER, otherwise you cannot run it by docker container.
1. git clone https://github.com/yuyuqi-design/IREDFisher.git
2. cd IREDFisher
3. cd data
4. #Download data from https://figshare.com/ndownloader/files/38466119 
   </br>
   wget -O IREDFisher_data.zip https://figshare.com/ndownloader/files/38466119
6. unzip IREDFisher_database.zip
7. cd ..
8. sudo docker-compose -f docker-compose.yml build --no-cache
9. sudo docker-compose -f docker-compose.yml up
10. open your favourite browser and paste http://127.0.0.1:5000 
11. login use email:admin@email.com and password as: password

### Enjoy!
