# Variant calling algorithm
You can run variant_calling.py scrypt with next command:
../variant_calling.py -v relative_path_to_pileup_file

Result of variant calling algortihm will be written in output/VCF_impl.vcf file.

In resources/VCF_ref.vcf file are referent results of variant calling executed with bcftools.

To compare results of implemented algortihm and bcftools, you need to run compare.py scrypt with next command:
../compare.py --compare relative_path_to_implemented_vcf relative_path_to_referent_vcf

Results of comparison are written in output/compare_results.txt file.

To run test, execute test.py scrypt with next command:
../test.py -t testType
where testType is number of test, with --help you can see list of tests.

All paths are relative in regard to folder which contains variant_calling.py scrypt.

Folder resources contains resources files (VCF, pileup).

For more details you can check presentation variant_calling_pm173351m.pptx

Also you can chek video presentation on youtube:
https://www.youtube.com/watch?v=kYzhIHWRGFk&feature=youtu.be
