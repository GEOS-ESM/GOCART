module ConvertUnits_mod
    use ConvertUnits_private_mod, only: read_scale_config, ScaleMapReal32, ScaleMapReal64

    implicit none
    private

    public read_scale_config
    public ScaleMapReal32
    public ScaleMapReal64
end module ConvertUnits_mod