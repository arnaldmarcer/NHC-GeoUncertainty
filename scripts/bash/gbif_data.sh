# ================================================================================================ #
# First download the GBIF derived dataset from:
# https://zenodo.org/record/5052596/files/20210311T1049_preserved_specimen_table.zip?download=1
# and placed it under a 'tmp' directory in your root directory.

# PROJECT ROOT DIRECTORY ------------------------------------------------------------------------- #
root_dir=<path>/NHC-GeoUncertainty # set your root directory and place zip file here
mkdir -p $root_dir/tmp
cd $root_dir/tmp

# UNZIP FILE
echo "Unzipping preserved specimens data ..."
unzip ../20210311T1049_preserved_specimen_table.zip

# SPLIT INTO CHUNKS TO PROCESS ------------------------------------------------------------------- #
echo "Splitting preserved specimens data into chunks of 1 000 000 records ..."
split -d -l 1000000 preserved_specimen_table_2.csv ps
cd $root_dir

# IMPORT INTO NEW SQLITE3 DATABASE --------------------------------------------------------------- #
db_filename="data/raw/gbif/20210311T1049_preserved_specimen_sqlite.db"
echo "Importing data into database "$db_filename" ..."
STARTTIME=$(date +%s)
for f in $root_dir/tmp/ps*;
do
    echo "  processing "$(basename $f)" ..."
    echo "    --> Striping quotes and backslashes from fields to prevent 'unescaped character error'..."
    cut -f 1-31 $f | sed 's/\"//g' | sed 's/\\//g' | sed 's/“//g' > $root_dir/tmp/tmp.csv
    echo "    --> Importing into $root_dir/$db_filename ..."
    sqlite3 -separator $'\t' $root_dir/$db_filename ".import tmp/tmp.csv ps"
    rm $root_dir/tmp/tmp.csv
    ENDTIME=$(date +%s)
done
# CREATE MAJOR INDEX TO SPLIT DATABASE WITH R SCRIPT --------------------------------------------- #
sqlite3 $root_dir/$db_filename "create index idx_hascoordinate on ps(hascoordinate);"
ENDTIME=$(date +%s)
echo "Seconds elapsed: $(($ENDTIME - $STARTTIME))"

# CLEAN TMP DIR
rm -r $root_dir/tmp
# ================================================================================================ #
