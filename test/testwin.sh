for f in obj/*.exe; do
  OUT=`wine $f > /dev/null` 
  if [ ! -z "$OUT" ]; then
    echo ERROR
    exit
  fi
done
echo OK!
