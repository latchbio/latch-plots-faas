#!/usr/bin/env bash
if [[ "$0" = "${BASH_SOURCE[0]}" ]];
then
  echo "Error: this script was not sourced."
  exit 1
fi

export domain='latch.bio'
export auto_reload='false'
export logging_mode='console'
