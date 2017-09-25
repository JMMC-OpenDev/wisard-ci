; catch any error to exit cleanly in batch mode (GDL does not - yet)
  if (wisard_is_interactive eq 0) then begin
     CATCH, error_status
     if (error_status ne 0) then begin
        CATCH,/CANCEL
        message,/reissue_last,/informational
        print,"non-interactive session: error has occured, exiting."
        exit, status=1
     end
  end
