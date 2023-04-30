/*
  Library that serves as the main interface for writing to an hdf5 dataset.
  The hierarchy employed is explain in /_readme/hdf5dataOverview.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

/* Macros */
#include "macros.h"
/* hdf5 library */
#include "hdf5.h"
/* File I/O functions */
#include "hdf5Interface.h" //main interface

/* Local function prototypes */
void hdf5InterfaceInitialise();
int touchFolderStructure(const char *fToTouch_ptr);
int countFolderNesting(const char *fToCount_ptr);
char** getFolderNesting(int folderCount, const char *fToCount_ptr);
int writeSimFolder(char *folderName);
int checkDatasetExists(const char *datToChk_ptr);
int checkFileExists(const char *fToChk_ptr);
hdf5Hierarchy getHdf5Data();
int writeDouble1DArr(double *dblArr_ptr,
		     int dblArr_sz,
		     const char *absolutePath,
		     const char *dataSetName);
double* readDouble1DArr(const char *absolutePath, const char *dataSetName);
/*
  Write a Complex 2D array to hdf5 file. [No chunking, compression etc.]
  Suppose we have A+IB, then data is written as two arrays:
  A with name: "dataname" + "-re"
  B with name: "dataname" + "-im"
  Absolute path is assumed to exist.
  If dataset exists, do not write.
  ret 1 if written, 0 if not.
*/
int writeComplexDouble2DArr(complex double **cplxdblArr_ptr,
			    int dblArr_szX, int dblArr_szY,
			    const char *absolutePath,
			    const char *dataSetName);
int writeUllint1DArr(ullint *ullArr_ptr,
		     ullint ullArr_sz,
		     const char *absolutePath,
		     const char *dataSetName);
ullint* readUllint1DArr(const char *absolutePath, const char *dataSetName);

//double* readDouble1DArr(const char *absolutePath, const char *dataSetName);


/*
  In this definition only the maximally nested object of interest need be
  defined.
*/
static const char *folderSystem[] = {"/simulationData", //ensure folder exists
				     "/simulationData/S2advection",
				     "/simulationData/S2advection/scalar",
				     "/simulationData/S2advection/vector",
				     "/calculationData", //ensure folder exists
				     "/calculationData/wigner",
				     "/calculationData/wigner/delFuncs",
				     "/calculationData/wigner/w3jFuncs"};

static const char fileName[] = "data.h5";
hdf5Hierarchy hdf5Info; //make info-struct global, init in main

/*
  To Do:
  Move (file)folderSystem to settings file
  complex 1d arrays
  Investigate hyperslabs/compression/Chunking
*/


/*
Populate structure with useful info on instantiation of library
 */
void hdf5InterfaceInitialise(){
  hdf5Info = getHdf5Data(); //populate info-struct
  touchFolderStructure(fileName); //ensure hdf5 file correctly initialised
}


/* 
   Ensure the target file exists, and establish the desired folder structure.
   -maybe add error checking.
*/
int touchFolderStructure(const char *fToTouch_ptr){
  
  hid_t file, grp;
  int i,j; //iterator(s)
 
  
  //Ensure that the file exists
  if(!checkFileExists(fToTouch_ptr)){
    file = H5Fcreate(fToTouch_ptr, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  file = H5Fopen(fToTouch_ptr, H5F_ACC_RDWR, H5P_DEFAULT);

  for(i=0; i<hdf5Info.folderSystem_sz; i++){
    //Total number of folders in current absolute path
    int fCou = countFolderNesting(folderSystem[i]);
    char **folders = getFolderNesting(fCou, folderSystem[i]);

    for(j=0; j<fCou; j++){
      //printf("%s, fc=%d\n", folders[j], fCou);
      
      if(H5Lexists(file, folders[j], H5P_DEFAULT)==0){
	//Group doesn't exist, create
	grp = H5Gcreate(file, folders[j], H5P_DEFAULT, H5P_DEFAULT, 
			H5P_DEFAULT);
      }
      
      free(folders[j]);
    }
    free(folders);
  }

  //Release
  H5Fclose(file);
  return 1;
}


/*
  Given a string "/folderA/folderB/folderC" count the number of folders
*/
int countFolderNesting(const char *fToCount_ptr){
  int i=0; //iterator(s)

  char *ch_ptr;
  char dirDelim ='/';

  ch_ptr=strchr(fToCount_ptr, dirDelim);

  while(ch_ptr!=NULL){ //return type of strchr is null_ptr if str not found
    //printf("folder[%d] delim @ %d\n", i,
    //   ch_ptr-fToCount_ptr+1);
    ch_ptr=strchr(ch_ptr+1, dirDelim);
    i++;
  }
 
  return i;
}


/*
  Given a string "/folderA/folderB/folderC" and a count of the number of folders
  return an array {"folderA", "folderA/folderB", "folderA/folderB/folderC"}
*/
char** getFolderNesting(int folderCount, const char *fToExtr_ptr){
  int i = 0; //iterator(s)

  char *ch_ptr;
  
  char dirDelim = '/';

  int fToExtr_sz = strlen(fToExtr_ptr);
  char strTmp[fToExtr_sz];

  char** fNames;
  fNames = (char**)malloc(folderCount*sizeof(char*));

  strncpy(strTmp, fToExtr_ptr+1, fToExtr_sz-1);
  strTmp[fToExtr_sz-1]= '/';

  ch_ptr=strchr(strTmp, dirDelim);
  while(ch_ptr!=NULL){ //return type of strchr is null_ptr if str not found
    
    int curLen = ch_ptr-strTmp;    
    
    if(i<folderCount){
      fNames[i] = (char*) malloc((1+curLen)*sizeof(char));
      
      strncpy(1+fNames[i], ch_ptr-curLen, curLen);
      fNames[i][0] = '/';
      //Must explicitly null-terminate when using strncpy!!
      fNames[i][curLen+1] = '\0';
      //printf("curStr=%s, %d, %d, s=%s\n", fNames[i], curLen,
      //       strlen(ch_ptr), ch_ptr);
      ch_ptr=strchr(ch_ptr+1, dirDelim);
    } else {
      ch_ptr=NULL;
    }
    
    i++;
  }

  return fNames;
}

/*
  Folder structure must already exist, here we are restricted to one 
  level of nesting deeper
*/
int writeSimFolder(char *folderName){  
  hid_t file, grp;
  file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
  grp = H5Gcreate(file, folderName, H5P_DEFAULT, H5P_DEFAULT, 
		  H5P_DEFAULT);
  //Release
  H5Fclose(file);
  return 1;
}

/*
  ret 1 if exists, 0 if not.
*/
int checkDatasetExists(const char *datToChk_ptr){
  hid_t file;
  file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

  if(H5Lexists(file, datToChk_ptr, H5P_DEFAULT)){
    //Group exists
    H5Fclose(file);
    return 1;
  }
  H5Fclose(file);
  return 0;
}



/* 
   ret 1 if exists, 0 if not
   Uses fopen, as this is portable. Could also use stat
*/
int checkFileExists(const char *fToChk_ptr){
  FILE *pFile_ptr;
  pFile_ptr = fopen(fToChk_ptr, "r");
  if(pFile_ptr==NULL){
    return 0; //file doesnt exist
  } else {
    fclose(pFile_ptr);
    return 1; //file exists
  }
}


/* Helper function to populate struct */
hdf5Hierarchy getHdf5Data(){
  hdf5Hierarchy h1;
  h1.folderSystem_ptr = folderSystem;
  h1.folderSystem_sz = numElem(folderSystem);
  return h1;
}


/*
  Write a 1D array to hdf5 file. [No chunking, compression etc.]
  Absolute path is assumed to exist.
  If dataset exists, do not write.
  ret 1 if written, 0 if not.
*/
int writeDouble1DArr(double *dblArr_ptr,
		     int dblArr_sz,
		     const char *absolutePath,
		     const char *dataSetName){

  hid_t file, dataset, dataspace, datatype;
  herr_t status;
  hsize_t arr_sz[1];

  //printf("sz:dblArr_ptr[dblArr_sz]=%f\n", dblArr_ptr[dblArr_sz-1]); 

  char strTmp[strlen(absolutePath)+1+strlen(dataSetName)];
  const char *dirDelim = "/";
  strcpy(strTmp, absolutePath);
  strcat(strTmp, dirDelim);
  strcat(strTmp, dataSetName);

  //check if dataset already exists
  if(!checkDatasetExists(strTmp)){

  
    //Open hdf5 file for R/W
    file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

  
    //rank=1, length of dataset, maximum length
    arr_sz[0] = dblArr_sz;
    dataspace = H5Screate_simple(1, arr_sz, NULL);

    //define datatype for the data in the file
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    //without compression
    dataset = H5Dcreate(file, strTmp, H5T_NATIVE_DOUBLE,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		      H5P_DEFAULT, dblArr_ptr);


    //close/release resources
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    
    H5Fclose(file);
    
    return 1; //dataset written
  }
  return 0; //dataset already exists
}

/*
  Read a 1f array from an hdf5 file. [No chunking, hyperslabs, compression etc.]
  Absolute path and dataset are assumed to exist.
*/
double* readDouble1DArr(const char *absolutePath, const char *dataSetName){

  hid_t       file, dataset;         /* handles */
  hid_t       datatype, dataspace;
  H5T_class_t t_class;                 /* data type class */
  hsize_t     dims_out[1];           /* dataset dimensions */
  herr_t      status;

  int status_n, rank;

  /* Construct dataset name from input */
  char strTmpDataName[strlen(absolutePath)+1+strlen(dataSetName)];
  const char *dirDelim = "/";
  strcpy(strTmpDataName, absolutePath);
  strcat(strTmpDataName, dirDelim);
  strcat(strTmpDataName, dataSetName);


  /*
   * Open the file and the dataset.
   */
  file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, strTmpDataName, H5P_DEFAULT);

  /*
   * Get datatype and dataspace handles and then query
   * dataset class, dimensions.
   */
  datatype  = H5Dget_type(dataset);    /* datatype handle */
  t_class     = H5Tget_class(datatype);

  /* Get Data Info */
  dataspace = H5Dget_space(dataset);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  
  double* retData = (double*) malloc(dims_out[0]*sizeof(double));
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		   retData);

  /*
   * Close/release resources.
   */
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Fclose(file);


  double *retData_ptr;
  retData_ptr = &retData[0];

  return retData_ptr;
}


/*
  Write a double 2D array to hdf5 file. [No chunking, compression etc.]
  Note that data must first be copied into a single contiguous block!
  Absolute path is assumed to exist.
  If dataset exists, do not write.
  ret 1 if written, 0 if not.
*/
int writeDouble2DArr(double **dblArr_ptr,
		     int dblArr_szX, int dblArr_szY,
		     const char *absolutePath,
		     const char *dataSetName){
  hid_t file, dataset, dataspace, datatype;
  herr_t status;
  hsize_t arr_sz[2];
  ullint i, j; //Iterator(s)

  char strTmp[strlen(absolutePath)+1+strlen(dataSetName)];

  const char *dirDelim = "/";
  sprintf(strTmp, "%s%s%s",
	  absolutePath, dirDelim, dataSetName);

  //check if dataset already exists
  if(!checkDatasetExists(strTmp)){

    //Open hdf5 file for R/W
    file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
  
    //rank=2, length of dataset, maximum length
    arr_sz[0] = dblArr_szX;
    arr_sz[1] = dblArr_szY;
    dataspace = H5Screate_simple(2, arr_sz, NULL);

    //define datatype for the data in the file
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    dataset = H5Dcreate(file, strTmp, H5T_NATIVE_DOUBLE,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //Prepare data for writing (must be contiguous)
    double** data = (double**) malloc(dblArr_szX*sizeof(double*));
    data[0] = (double*) malloc(dblArr_szX*dblArr_szY*sizeof(double));
    
    //Update pointers
    for(i=1; i<dblArr_szX; i++){
      data[i] = data[0]+i*dblArr_szY;
    }

    //Populate arrays
    for(i=0; i<dblArr_szX; i++){
      for(j=0; j<dblArr_szY; j++){
	data[i][j]=creal(dblArr_ptr[i][j]);
      }
    }

    //Write
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		      H5P_DEFAULT, &data[0][0]);

    //close/release hdf resources
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    
    H5Fclose(file);
    
    //Free copied data
    free(data[0]);
    free(data);
    return 1; //dataset written
  }

  return 0;

}


/*
  Write a Complex 2D array to hdf5 file. [No chunking, compression etc.]
  Suppose we have A+IB, then data is written as two arrays:
  A with name: "dataname" + "-re"
  B with name: "dataname" + "-im"
  Note that data must first be copied into a single contiguous block!
  Absolute path is assumed to exist.
  If dataset exists, do not write.
  ret 1 if written, 0 if not.
*/
int writeComplexDouble2DArr(complex double **cplxdblArr_ptr,
			    int dblArr_szX, int dblArr_szY,
			    const char *absolutePath,
			    const char *dataSetName){
  
  hid_t file, datasetRe, datasetIm, dataspace, datatype;
  herr_t statusRe, statusIm;
  hsize_t arr_sz[2];
  ullint i, j; //Iterator(s)

  char strTmpR[strlen(absolutePath)+1+strlen(dataSetName)+3];
  char strTmpI[strlen(absolutePath)+1+strlen(dataSetName)+3];
  
  const char *dirDelim = "/";
  sprintf(strTmpR, "%s%s%s-Re",
	  absolutePath, dirDelim, dataSetName);
  sprintf(strTmpI, "%s%s%s-Im",
	  absolutePath, dirDelim, dataSetName);

  //check if dataset already exists
  if((!checkDatasetExists(strTmpR)) && (!checkDatasetExists(strTmpI))){

    //Open hdf5 file for R/W
    file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

  
    //rank=2, length of dataset, maximum length
    arr_sz[0] = dblArr_szX;
    arr_sz[1] = dblArr_szY;
    dataspace = H5Screate_simple(2, arr_sz, NULL);

    //define datatype for the data in the file
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    statusRe = H5Tset_order(datatype, H5T_ORDER_LE);
    statusIm = H5Tset_order(datatype, H5T_ORDER_LE);

    //Real part of data (without compression)
    datasetRe = H5Dcreate(file, strTmpR, H5T_NATIVE_DOUBLE,
			  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    //Imag part of data (without compression)
    datasetIm = H5Dcreate(file, strTmpI, H5T_NATIVE_DOUBLE,
			  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    //Prepare data for writing (must be contiguous)
    double** dataRe = (double**) malloc(dblArr_szX*sizeof(double*));
    dataRe[0] = (double*) malloc(dblArr_szX*dblArr_szY*sizeof(double));
    
    double** dataIm = (double**) malloc(dblArr_szX*sizeof(double*));
    dataIm[0] = (double*) malloc(dblArr_szX*dblArr_szY*sizeof(double));

    //Update pointers
    for(i=1; i<dblArr_szX; i++){
      dataRe[i] = dataRe[0]+i*dblArr_szY;
      dataIm[i] = dataIm[0]+i*dblArr_szY;
    }

    //Populate arrays
    for(i=0; i<dblArr_szX; i++){
      for(j=0; j<dblArr_szY; j++){
	dataRe[i][j]=creal(cplxdblArr_ptr[i][j]);
	dataIm[i][j]=cimag(cplxdblArr_ptr[i][j]);
      }
    }

    //Write
    statusRe = H5Dwrite(datasetRe, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		      H5P_DEFAULT, &dataRe[0][0]);
    statusIm = H5Dwrite(datasetIm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		      H5P_DEFAULT, &dataIm[0][0]);

    //close/release hdf resources
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(datasetRe);
    H5Dclose(datasetIm);

    H5Fclose(file);
    
    //Kill copied arrs
    //Free real portion of copied data
    free(dataRe[0]);
    free(dataRe);
    //Free imag portion of copied data
    free(dataIm[0]);
    free(dataIm);

    return 1; //dataset written
  }

  return 0;
}


/*
  Write a 1D array to hdf5 file. [No chunking, compression etc.]
  Absolute path is assumed to exist.
  If dataset exists, do not write.
  ret 1 if written, 0 if not.
*/
int writeUllint1DArr(ullint *ullArr_ptr,
		     ullint ullArr_sz,
		     const char *absolutePath,
		     const char *dataSetName){

  hid_t file, dataset, dataspace, datatype;
  herr_t status;
  hsize_t arr_sz[1];

  //printf("sz:ullArr_ptr[ullArr_sz]=%llu\n", ullArr_ptr[ullArr_sz-1]); 

  char strTmp[strlen(absolutePath)+1+strlen(dataSetName)];
  const char *dirDelim = "/";
  strcpy(strTmp, absolutePath);
  strcat(strTmp, dirDelim);
  strcat(strTmp, dataSetName);

  //check if dataset already exists
  if(!checkDatasetExists(strTmp)){

  
    //Open hdf5 file for R/W
    file = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

  
    //rank=1, length of dataset, maximum length
    arr_sz[0] = ullArr_sz;
    dataspace = H5Screate_simple(1, arr_sz, NULL);

    //define datatype for the data in the file
    datatype = H5Tcopy(H5T_NATIVE_ULLONG);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    //without compression
    dataset = H5Dcreate(file, strTmp, H5T_NATIVE_ULLONG,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, 
		      H5P_DEFAULT, ullArr_ptr);


    //close/release resources
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    
    H5Fclose(file);
    
    return 1; //dataset written
  }
  return 0; //dataset already exists
}

/*
  Read a ullint array from an hdf5 file.
  [No chunking, hyperslabs, compression etc.]
  Absolute path and dataset are assumed to exist.
*/
ullint* readUllint1DArr(const char *absolutePath, const char *dataSetName){

  hid_t       file, dataset;         /* handles */
  hid_t       datatype, dataspace;
  H5T_class_t t_class;                 /* data type class */
  hsize_t     dims_out[1];           /* dataset dimensions */
  herr_t      status;

  int status_n, rank;

  /* Construct dataset name from input */
  char strTmpDataName[strlen(absolutePath)+1+strlen(dataSetName)];
  const char *dirDelim = "/";
  strcpy(strTmpDataName, absolutePath);
  strcat(strTmpDataName, dirDelim);
  strcat(strTmpDataName, dataSetName);


  /*
   * Open the file and the dataset.
   */
  file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(file, strTmpDataName, H5P_DEFAULT);

  /*
   * Get datatype and dataspace handles and then query
   * dataset class, dimensions.
   */
  datatype  = H5Dget_type(dataset);    /* datatype handle */
  t_class     = H5Tget_class(datatype);

  /* Get Data Info */
  dataspace = H5Dget_space(dataset);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  
  ullint* retData = (ullint*) malloc(dims_out[0]*sizeof(ullint));
  status = H5Dread(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		   retData);

  /*
   * Close/release resources.
   */
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Fclose(file);


  ullint *retData_ptr;
  retData_ptr = &retData[0];

  return retData_ptr;
}
