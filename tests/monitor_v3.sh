#!/bin/bash
echo "#=================================="
echo "#=========== MONITORING ==========="
echo "#=================================="
echo "#--- General Information ---"
echo "#CPU: $(nproc)"
echo "#Total Memory: $(free -h | grep Mem | awk '{ print $2 }')"
echo "#Total Disk space: $(df -h | grep cromwell_root | awk '{ print $2}')"
echo "#"


echo "#time cpu_usage_percent mem_avail_gb disk_free_gb"

function runtimeInfo() {
        datetime=$(date "+%F_%T")
        cpu=$(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')
        mem_avail=$(free -g | grep Mem | awk '{print $7}')
        disk_free=$(df -BG | grep cromwell | awk '{ print $4 }' | tr -d G)
        echo $datetime $cpu $mem_avail $disk_free
}

while true; do runtimeInfo; sleep 60; done
