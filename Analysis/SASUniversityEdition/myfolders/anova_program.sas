*initialize empty results table;
proc sql;
	create table alldata  (MS_sp num, MS_st num);
	run;
*macro that performs the analysis given gene name and total results table;
%macro modelrandom(gene_name);

FILENAME REFFILE """/folders/myfolders/test_input_data_male/&gene_name\.txt""";
proc import datafile=REFFILE
	out = sample_data
	DBMS=dlm
	replace;
	run;

ods output Type3=anova_table;
proc mixed data=sample_data cl method=type3;

     class species strains;
     model exp= ;
     random species strains(species);
     run;

data ms_table;
	set anova_table;
	keep Source MS;
	if Source = 'Residual' then delete;
run;
proc transpose data=ms_table
	out = ms_table_transposed;
run;
data ms_table_transposed;
	set ms_table_transposed;
	keep COL1 COL2;
	rename COL1=MS_sp;
	rename COL2=MS_st;
	
run;


proc append base=alldata data=ms_table_transposed ;
	run;
%mend;




proc import datafile="/folders/myfolders/test_input_data_male/genelist.txt"
	out = genelist
	DBMS=dlm
	replace;
	run;
	

*run macros over entire gene list;	
data run_macros;
	set genelist;
	str = catt('%modelrandom(',gene,');');
	call execute (str);

run;

*Save the alldata output;

proc export data=alldata 
      outfile="/folders/myfolders/anova_output/test_anova_output_male.txt"
      dbms=dlm;  
      delimiter=' ';
run;

proc export data=genelist 
      outfile="/folders/myfolders/anova_output/test_genelist_male.txt"
      dbms=dlm;  
      delimiter=' ';
run;
	

