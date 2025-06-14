for i in {1..100}; do
  echo "Deleting CSVs in output_$i …"
  find "./output_${i}" -maxdepth 1 -type f -name '*.csv' \
       -newermt '2025-04-12' ! -newermt '2025-04-14' \
       -delete
done
