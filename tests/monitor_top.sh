#!/bin/bash
echo "#==========================================#"
echo "#=========== PROCESS MONITORING ===========#"
echo "#==========================================#"

echo "#time cpu_usage_percent mem_usage_percent process"

function runtimeInfo() {
        datetime=$(date "+%F_%T")
		top -b -n1 | grep -A100 COMMAND | grep -v COMMAND | awk -v dt=$datetime 'length>20 { print dt,$9,$10,$12 }'
}

while true; do runtimeInfo; sleep 60; done