# Variant calling algorithm
You need to run variant_calling.py scrypt with appropriate option:
-v is used for variant calling algorithm,
-t is used for tests, and you also need to enter number of test,
-c is used for comparison and
-h is for help.

Usage: 

       variant_calling.py [-vch]  ...
       
       variant_calling.py [-th] [testType] ...
       
testType:

  1 variantCallingCore test
  
  2 variantCalling test
  
  3 addVarinatToIntVcArray test
  
Options:

  -h       help 
  
  -v       execute variant calling
  
  -t Argument       execute tests
  
  -c       compare results

Result of variant calling algortihm will be written in output/VCF_impl.vcf file.

In resources/VCF_ref.vcf file are referent results of variant calling executed with bcftools.

To compare results of implemented algortihm and bcftools, you need to run compare option, 
it uses output/VCF_impl.vcf and resources/VCF_ref.vcf files, and they need to be in the same folder as this scrypt.

Results of comparison are written in output/compare_results.txt file.

For more details you can check presentation variant_calling_pm173351m.pptx

Also you can chek video presentation on youtube:
https://www.youtube.com/watch?v=kYzhIHWRGFk&feature=youtu.be
