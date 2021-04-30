#!/bin/bash

record_bin="glc-capture"
real_graphics_bin="milkyway_nbody_graphics"
ffmpeg_bin="ffmpeg"

command -v ${record_bin} >/dev/null 2>&1 || { echo >&2 "${record_bin} not found"; exit 1; }
command -v ${real_graphics_bin} >/dev/null 2>&1 || { echo >&2 "${real_graphics_bin} not found"; exit 1; }
command -v ${real_graphics_bin} >/dev/null 2>&1 || { echo >&2 "${ffmpeg_bin} not found"; exit 1; }


width=1024
height=768

output_name="nbody_video"
output_file="${output_name}.glc"
video_codec="libx264"
output_video_file="${output_name}.mp4"

if [ -e ${output_video_file} ]; then
    echo "Output video '${output_video_file}' already exists"
    exit 1;
fi

echo "Launching graphics. Recording to ${output_file}"

${record_bin} --start                                   \
              --fps=30                                  \
              --disable-audio                           \
              --out ${output_file}                      \
               ${real_graphics_bin} $*                  \
                                    --width=${width}    \
                                    --height=${height}  \
                                    --block-simulation  \
                                    --quit-on-complete  \
                                    --no-show-info

echo "Recording done. Beginning encoding video."

glc-play ${output_file} -o - -y 1 | \
    ffmpeg -i - -an -vcodec ${video_codec} -crf 22 -threads 0 ${output_video_file}


if [ $? -eq 0 ]; then
    echo "Video encoding done: ${output_video_file}"
else
    echo "Error encoding video"
    exit 1
fi

