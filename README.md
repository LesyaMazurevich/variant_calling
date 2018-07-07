# Variant calling algorithm

After running variant_calling.py scrypt, you need to enter appropriate value to continue:

1 - execute variant calling algorithm

2 - execute variantCallingCore function test

3 - execute variantCalling function test

4 - execute addVariantToIncArray function test

Result of variant calling algortihm will be written in VCF_impl.vcf file.

In VCF_ref.vcf file are referent results of variant calling executed with bcftools.

To compare results of implemented algortihm and bcftools, you need to run compare.py scrypt, 
it uses VCF_impl.vcf and VCF_ref.vcf files, and they need to be in the same folder as this scrypt.

Results of comparison are written in compare_results.txt file.

For more details you can check presentation variant_calling_pm173351m.pptx

Also you can chek video presentation on youtube:
https://www.youtube.com/watch?v=kYzhIHWRGFk&feature=youtu.be
