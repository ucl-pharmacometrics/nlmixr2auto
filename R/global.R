# ================================================================
# Global variable declarations
# Each section corresponds to a specific function module.
# ================================================================

if (getRversion() >= "3.5.0") {
  # ------------------------ mod.run ----------------------------
  utils::globalVariables(c(
    "f"
  ))
  # ------------------------ modelGen ----------------------------
  utils::globalVariables(c("Name", "init"))

  # ------------------------ sf.operator ----------------------------
  utils::globalVariables(c("round.num"))

}
