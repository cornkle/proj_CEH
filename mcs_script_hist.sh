#!/bin/bash

for i in {16..23}
do
  ipython GLOBAL/saveMCS_CP4_storm.py $i CP4hist
done
