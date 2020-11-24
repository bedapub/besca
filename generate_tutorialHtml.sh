# GET static file:
for file in adata_to_eset.ipynb bescape_tutorial.ipynb
do
    jupyter nbconvert --to html ../tutorials/$file 
done
mv ../tutorials/*.html .
