#!/bin/sh

export NVM_DIR="$HOME/.nvm"
source "$NVM_DIR/nvm.sh"
nvm install --lts
nvm use --lts
npm install -g @anthropic-ai/claude-code
HOME=`pwd`
