#!/usr/bin/env bash

# Shared logging helpers for Nextflow module scripts.

log_info() {
  echo "[INFO]: ℹ️ $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warning() {
  echo "[WARNING]: ⚠️ $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_error() {
  echo "[ERROR]: ❌ $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
  exit 1
}

log_success() {
  echo "[SUCCESS]: ✅✅✅ $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_save() {
    echo "[SAVE]: 💾 $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_start() {
    echo "[START]: 🚀 $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_ok() {
    echo "[OK]: 🟢 $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_verify() {
    echo "[VERIFY]: 🔍 $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_log() {
    echo "[LOG]: 📝 [$(date '+%Y-%m-%d %H:%M:%S')] $1"
}