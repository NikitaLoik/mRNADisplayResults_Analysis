# 1 GET NAMES =========================================================

# 1.1 Get virtual-environment name.
VENV_NAME := $(wildcard .venv_*)

# 1.2 Get virtual-environment stem.
VENV_STEM := $(subst .venv_,,$(VENV_NAME))

# 1.3 Get kernel names from virtual-environment stem.
KERNEL_NAME_PY = "$(VENV_STEM)_py"

# =====================================================================

venv_get:
	scripts_shell/get_venv.sh

venv_remove:
	@echo "===removing virtual environment==="
	rm -rf $(VENV_NAME)

	@echo "===removing virtual iPython kernel==="
	rm -rf ~/Library/Jupyter/kernels/$(KERNEL_NAME_PY)
