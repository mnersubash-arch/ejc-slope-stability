Attribute VB_Name = "Thesis_Columns"
Option Explicit
' GENG5512 thesis package v14.
' Primary VBA Hungr simplified column assembly and genuine central-section 2D Bishop assembly.
' Source implementation supplied for native import/compile/run. NOT executed in Microsoft Excel here.
' Reads Interface inputs; builds its OWN variable-resolution grid; never reads worksheet strength/F values.
' Strength remains eq. 4.21, with analytical Newton on phi via EJC_Eq421. No exact-EJC C1 substitution.
' Closure 3 is explicitly a conventional-HB comparator. MC is an explicit separate verification model.
' Geometry is restricted to the documented positive-centre, driving-side, toe-cutoff model.
Private Const PI_C As Double = 3.141592653589793
Private Type ModelConfig
    xc As Double
    zc As Double
    radius As Double
    height As Double
    slope As Double
    gamma As Double
    uExtra As Double
    waterAngle As Double
    waterHeight As Double
    strength As Long
    closure As Long
    nu As Double
    sci As Double
    mb As Double
    s As Double
    a As Double
    cMC As Double
    phiMC As Double
    tolA As Double
    tolE As Double
    tolF As Double
    checks As Long
    maxIter As Long
End Type
Private Type ModelElement
    x As Double
    y As Double
    zb As Double
    h As Double
    area As Double
    ah As Double
    g As Double
    kap As Double
    rho As Double
    weight As Double
    pressure As Double
    sn As Double
    phi As Double
    cohesion As Double
    finalN As Double
    tau As Double
    angular As Double
    ejc As Double
    closureRes As Double
    equilibrium As Double
    s1 As Double
    s2 As Double
    s3 As Double
    label As String
End Type
Private Type ModelResult
    status As String
    reason As String
    fos As Double
    moment As Double
    iterations As Long
    active As Long
    admissible As Long
    evaluated As Boolean
    dF As Double
    globalRes As Double
    maxAngle As Double
    maxEJC As Double
    maxEq As Double
    minNormal As Double
    tension As Long
    projections As Long
    volume As Double
End Type

Private Function NumberAt(ByVal ws As Worksheet, ByVal addr As String) As Double
    Dim v As Variant
    v = ws.Range(addr).Value2
    If IsError(v) Or IsEmpty(v) Then Err.Raise 5, , "Missing/invalid numeric input: " & addr
    If VarType(v) = vbBoolean Or VarType(v) = vbString Then Err.Raise 5, , "Numeric input required: " & addr
    If Not IsNumeric(v) Then Err.Raise 5, , "Numeric input required: " & addr
    NumberAt = CDbl(v)
End Function

Private Function IntegerAt(ByVal ws As Worksheet, ByVal addr As String) As Long
    Dim v As Double
    v = NumberAt(ws, addr)
    If v <> Fix(v) Then Err.Raise 5, , "Integer input required: " & addr
    IntegerAt = CLng(v)
End Function

Private Function GeometryIssue(ByRef cfg As ModelConfig) As String
    Dim xToe As Double
    If cfg.radius <= cfg.zc Or cfg.zc <= 0 Or cfg.height <= 0 Or cfg.slope <= 0 Or cfg.gamma <= 0 Then
        GeometryIssue = "UNSUPPORTED GEOMETRY: require R > zc > 0 and positive height, slope and unit weight."
        Exit Function
    End If
    If cfg.xc <= 0 Or cfg.xc >= cfg.radius Then
        GeometryIssue = "UNSUPPORTED GEOMETRY: this driving-side toe-cutoff model requires 0 < xc < R."
        Exit Function
    End If
    xToe = cfg.xc - Sqr(cfg.radius ^ 2 - cfg.zc ^ 2)
    ' Permit the published reference circle's 2.62 mm rounding mismatch; do not move its grid.
    If xToe < -0.00001 * cfg.radius Then
        GeometryIssue = "UNSUPPORTED GEOMETRY: the adopted grid would omit slope-side mass between its right bound and the toe."
    End If
End Function

Private Sub ReadConfig(ByRef cfg As ModelConfig)
    Dim ws As Worksheet, issue As String
    Set ws = ThisWorkbook.Worksheets("Interface")
    If NumberAt(ws, "D6") <> 0 Then Err.Raise 5, , "Sphere centre y must be zero for the symmetry model."
    cfg.xc = NumberAt(ws, "C6"): cfg.zc = NumberAt(ws, "E6"): cfg.radius = NumberAt(ws, "C9")
    cfg.height = NumberAt(ws, "D14"): cfg.slope = NumberAt(ws, "D16") / NumberAt(ws, "E16")
    cfg.gamma = NumberAt(ws, "D22"): cfg.uExtra = NumberAt(ws, "D23")
    cfg.waterAngle = NumberAt(ws, "D24"): cfg.waterHeight = NumberAt(ws, "D25")
    cfg.strength = IntegerAt(ws, "P4"): cfg.closure = IntegerAt(ws, "P5"): cfg.nu = NumberAt(ws, "P6")
    cfg.sci = NumberAt(ws, "P7")
    cfg.mb = NumberAt(ws, "P8") * Exp((NumberAt(ws, "P9") - 100) / (28 - 14 * NumberAt(ws, "P10")))
    cfg.s = Exp((NumberAt(ws, "P9") - 100) / (9 - 3 * NumberAt(ws, "P10")))
    If IsError(ws.Range("P14").Value2) Then Err.Raise 5, , "Invalid a input."
    If IsNumeric(ws.Range("P14").Value2) And VarType(ws.Range("P14").Value2) <> vbBoolean Then
        cfg.a = CDbl(ws.Range("P14").Value2)
    Else
        cfg.a = 0.5 + (Exp(-NumberAt(ws, "P9") / 15) - Exp(-20# / 3)) / 6
    End If
    cfg.cMC = NumberAt(ws, "D20"): cfg.phiMC = NumberAt(ws, "D21")
    cfg.tolA = NumberAt(ws, "P16"): cfg.tolE = NumberAt(ws, "P17"): cfg.tolF = NumberAt(ws, "P18")
    cfg.checks = IntegerAt(ws, "P19"): cfg.maxIter = 2000
    If cfg.checks <> 0 And cfg.checks <> 1 Then Err.Raise 5, , "Acceptance switch must be 0 or 1."
    issue = GeometryIssue(cfg): If Len(issue) > 0 Then Err.Raise 5, , issue
    If cfg.strength < 0 Or cfg.strength > 1 Or cfg.closure < 1 Or cfg.closure > 3 Then Err.Raise 5, , "Unsupported model/closure."
    If cfg.sci <= 0 Or cfg.mb <= 0 Or cfg.s < 0 Or cfg.a <= 0 Or cfg.a > 1 Or cfg.nu <= 0 Or cfg.nu >= 0.5 Then Err.Raise 5, , "Invalid rock/Poisson input."
    If cfg.cMC < 0 Or cfg.phiMC < 0 Or cfg.phiMC >= 90 Or cfg.uExtra < 0 Or cfg.waterHeight < 0 Or cfg.waterAngle < 0 Or cfg.waterAngle >= 90 Then Err.Raise 5, , "Invalid strength/water input."
    If cfg.tolA <= 0 Or cfg.tolA > 0.0000000001 Or cfg.tolE <= 0 Or cfg.tolE > 0.00000001 Or cfg.tolF <= 0 Or cfg.tolF > 0.00000001 Then Err.Raise 5, , "Invalid acceptance tolerances."
End Sub

Private Sub ClosureCoefficients(ByRef cfg As ModelConfig, ByRef pC As Double, ByRef qC As Double)
    Select Case cfg.closure
        Case 1: pC = -1: qC = 1 / cfg.nu
        Case 2: pC = cfg.nu: qC = cfg.nu
        Case 3: pC = 0: qC = 1
    End Select
End Sub

Private Sub BuildElements(ByRef cfg As ModelConfig, ByVal dimN As Long, ByVal meshN As Long, ByRef els() As ModelElement, ByRef n As Long, ByRef volume As Double)
    Dim i As Long, j As Long, ny As Long, x0 As Double, x1 As Double, dx As Double, dy As Double
    Dim xx As Double, yy As Double, q As Double, base As Double, top As Double, thick As Double, water As Double, ah As Double
    Dim issue As String
    issue = GeometryIssue(cfg): If Len(issue) > 0 Then Err.Raise 5, , issue
    If (dimN <> 2 And dimN <> 3) Or meshN < 2 Or meshN > 20000 Then Err.Raise 5, , "Invalid mesh."
    If dimN = 3 And meshN > 400 Then Err.Raise 5, , "3D limit is 400 x 400; reduce mesh."
    ny = 1: If dimN = 3 Then ny = meshN
    ReDim els(1 To meshN * ny)
    x0 = cfg.xc - Sqr(cfg.radius ^ 2 - cfg.zc ^ 2): x1 = cfg.xc - cfg.radius
    dx = (x0 - x1) / meshN: dy = 1
    If dimN = 3 Then dy = cfg.radius / meshN
    n = 0: volume = 0
    For i = 1 To meshN
        xx = x0 - (i - 0.5) * dx
        For j = 1 To ny
            yy = 0: If dimN = 3 Then yy = (j - 0.5) * dy
            q = cfg.radius ^ 2 - (xx - cfg.xc) ^ 2 - yy ^ 2
            If q > 0 And xx <= 0 Then
                base = cfg.zc - Sqr(q): top = -xx * cfg.slope
                If top > cfg.height Then top = cfg.height
                thick = top - base
                If thick > 0 Then
                    n = n + 1: ah = dx * dy
                    With els(n)
                        .x = xx: .y = yy: .zb = base: .h = thick: .ah = ah
                        .g = Sqr(q) / cfg.radius: .area = ah / .g
                        .rho = Sqr(cfg.radius ^ 2 - yy ^ 2): .kap = (cfg.xc - xx) / .rho
                        .weight = cfg.gamma * thick * ah: .pressure = cfg.uExtra
                        If cfg.waterHeight > 0 Then
                            water = -xx * Tan(cfg.waterAngle * PI_C / 180)
                            If water > cfg.waterHeight Then water = cfg.waterHeight
                            If water > base Then .pressure = .pressure + 9.81 * (water - base)
                        End If
                        .sn = (cfg.gamma * thick - .pressure) * .g ^ 2: .phi = -1
                        .label = "NOT EVALUATED"
                    End With
                    volume = volume + thick * ah
                End If
            End If
        Next j
    Next i
    If n = 0 Then Err.Raise 5, , "No active elements."
    ReDim Preserve els(1 To n)
    If dimN = 3 Then volume = 2 * volume
End Sub

Private Function SolveStrength(ByRef cfg As ModelConfig, ByRef el As ModelElement, ByVal pC As Double, ByVal qC As Double) As Boolean
    Dim ph As Variant, phi As Double, hT As Double, b As Double, solverTol As Double
    If cfg.strength = 0 Then
        el.phi = cfg.phiMC: el.cohesion = cfg.cMC: SolveStrength = True: Exit Function
    End If
    solverTol = cfg.tolA
    If pC = 0 And qC = 1 Then solverTol = Application.Min(0.000000000001, cfg.tolA / 10)
    ph = EJC_FASTPHI(el.sn, cfg.sci, cfg.mb, cfg.s, cfg.a, pC, qC, el.phi, solverTol)
    If IsError(ph) Then el.label = "NO CONVERGED EQ.4.21 ROOT": Exit Function
    phi = CDbl(ph) * PI_C / 180: b = cfg.mb * el.sn / (1000 * cfg.sci) + cfg.s
    If b <= 0 Then el.label = "OUTSIDE HB DOMAIN": Exit Function
    hT = (b / (1 + Sin(phi) / cfg.a)) ^ cfg.a
    el.phi = CDbl(ph)
    el.cohesion = 1000 * cfg.sci * hT * Cos(phi) / 2 - el.sn * Tan(phi)
    el.label = EJC_ACCEPT_STATE(el.sn, el.cohesion, el.phi, cfg.sci, cfg.mb, cfg.s, cfg.a, pC, qC, cfg.tolA, cfg.tolE)
    SolveStrength = (el.label = "ADMISSIBLE")
End Function

Private Function ValidateFinal(ByRef cfg As ModelConfig, ByRef els() As ModelElement, ByVal n As Long, ByVal fos As Double, ByRef result As ModelResult, Optional ByVal rebuild As Boolean = True) As Boolean
    Dim i As Long, t As Double, ma As Double, normal As Double, label As String, chk As Variant, pC As Double, qC As Double
    Dim resist As Double, drive As Double, momentR As Double, momentD As Double
    ClosureCoefficients cfg, pC, qC
    result.admissible = 0: result.evaluated = True: result.maxAngle = 0: result.maxEJC = 0: result.maxEq = 0
    result.tension = 0: result.minNormal = 1E+300
    For i = 1 To n
        With els(i)
            t = Tan(.phi * PI_C / 180): ma = .g + .kap * t / fos
            label = "ADMISSIBLE"
            If ma <= 0 Or .g <= 0 Or .weight <= 0 Then
                .label = "INVALID EQUILIBRIUM": result.evaluated = False: ValidateFinal = False: Exit Function
            End If
            If rebuild Then .finalN = (.weight - .pressure * .ah - .cohesion * .area * .kap / fos) / (.area * ma)
            normal = .finalN: .tau = .cohesion + normal * t
            .equilibrium = Abs(.weight - ((normal + .pressure) * .ah + .tau * .area * .kap / fos)) / .weight
            If .equilibrium > result.maxEq Then result.maxEq = .equilibrium
            If normal < result.minNormal Then result.minNormal = normal
            If normal < 0 Then result.tension = result.tension + 1
            If .tau < 0 Or .equilibrium > 0.000000001 Then label = "REJECTED: EQUILIBRIUM/SHEAR"
            If cfg.strength = 1 Then
                chk = EJC_CHECK_STATE(normal, .cohesion, .phi, cfg.sci, cfg.mb, cfg.s, cfg.a, pC, qC, cfg.tolA)
                .label = EJC_ACCEPT_STATE(normal, .cohesion, .phi, cfg.sci, cfg.mb, cfg.s, cfg.a, pC, qC, cfg.tolA, cfg.tolE)
                If .label <> "ADMISSIBLE" Then label = .label
                If chk(7) = 1 Then
                    .angular = chk(1): .ejc = chk(2): .closureRes = chk(3): .s1 = chk(5): .s3 = chk(6): .s2 = pC * .s1 + qC * .s3
                    If Abs(.angular) > result.maxAngle Then result.maxAngle = Abs(.angular)
                    If Abs(.ejc) > result.maxEJC Then result.maxEJC = Abs(.ejc)
                Else
                    result.evaluated = False
                End If
            End If
            .label = label
            If label = "ADMISSIBLE" Then result.admissible = result.admissible + 1
            resist = resist + .tau * .area: drive = drive + .weight * .kap
            momentR = momentR + .tau * .area * .rho: momentD = momentD + .weight * .kap * .rho
        End With
    Next i
    If drive <= 0 Or momentD <= 0 Then ValidateFinal = False: Exit Function
    result.globalRes = Abs(resist / drive - fos) / fos: result.moment = momentR / momentD
    ValidateFinal = (result.evaluated And result.admissible = n And result.globalRes <= cfg.tolF)
End Function

Private Sub SolveModel(ByRef cfg As ModelConfig, ByVal dimN As Long, ByVal meshN As Long, ByRef els() As ModelElement, ByRef result As ModelResult, Optional ByVal logOutput As Boolean = False)
    Dim n As Long, i As Long, it As Long, pC As Double, qC As Double, fos As Double, fn As Double
    Dim t As Double, ma As Double, sf As Double, resist As Double, drive As Double, raw As Double, maxMove As Double
    Dim limit As Double, failures As Long, minA As Double, issue As String
    result = EmptyResult()
    result.status = "NOT ACCEPTED": result.reason = "NOT RUN"
    issue = GeometryIssue(cfg)
    If Len(issue) > 0 Then result.reason = issue: Erase els: Exit Sub
    BuildElements cfg, dimN, meshN, els, n, result.volume: result.active = n
    ClosureCoefficients cfg, pC, qC
    limit = -1000 * cfg.s * cfg.sci / cfg.mb: fos = 1
    If cfg.strength = 1 Then
        For i = 1 To n
            If els(i).sn <= limit Then els(i).sn = limit + Application.Max(0.000000001, 0.001 * Abs(limit)): result.projections = result.projections + 1
        Next i
    End If
    For it = 1 To cfg.maxIter
        resist = 0: drive = 0: failures = 0
        For i = 1 To n
            If Not SolveStrength(cfg, els(i), pC, qC) Then failures = failures + 1
            drive = drive + els(i).weight * els(i).kap
        Next i
        If failures > 0 Then
            result.iterations = it: result.reason = CStr(failures) & " of " & n & " elements have no admissible eq.4.21 strength; no FoS."
            Exit Sub
        End If
        If drive <= 0 Then result.reason = "Nonpositive driving sum.": Exit Sub
        maxMove = 0
        For i = 1 To n
            With els(i)
                t = Tan(.phi * PI_C / 180): ma = .g + .kap * t / fos
                If ma <= 0 Then result.reason = "Nonpositive m_alpha.": Exit Sub
                sf = (.cohesion * .ah + (.weight - .pressure * .ah) * t) / ma
                resist = resist + sf
                raw = (.weight - .pressure * .ah - .cohesion * .area * .kap / fos) / (.area * ma)
                If cfg.strength = 1 And raw <= limit Then
                    raw = limit + Application.Max(0.000000001, 0.001 * Abs(limit)): result.projections = result.projections + 1
                End If
                If Abs(raw - .sn) > maxMove Then maxMove = Abs(raw - .sn)
                .sn = raw
            End With
        Next i
        fn = resist / drive
        If fn <= 0 Then result.reason = "Nonpositive trial FoS.": Exit Sub
        result.dF = Abs(fn - fos) / fn: result.iterations = it: result.fos = fn
        If logOutput Then WriteIteration dimN, it, fn, result.dF, maxMove, result.projections
        If result.dF <= Application.Min(0.000000000001, cfg.tolF / 10) And maxMove <= 0.000000001 Then
            If ValidateFinal(cfg, els, n, fn, result) Then
                result.status = "ACCEPTED": result.reason = "Every final element and global equilibrium passed."
                If cfg.checks <> 1 Then result.status = "DIAGNOSTIC ONLY": result.reason = "Acceptance switch is off; displayed FoS withheld."
                Exit Sub
            End If
        End If
        fos = fn
        If it Mod 10 = 0 Then
            Application.StatusBar = "VBA " & dimN & "D: iteration " & it & ", mesh " & meshN
            DoEvents
        End If
    Next it
    result.reason = "Iteration limit reached; no accepted FoS."
End Sub

Private Function GetSheet(ByVal nameIn As String) As Worksheet
    On Error Resume Next
    Set GetSheet = ThisWorkbook.Worksheets(nameIn)
    On Error GoTo 0
    If GetSheet Is Nothing Then
        Set GetSheet = ThisWorkbook.Worksheets.Add(After:=ThisWorkbook.Worksheets(ThisWorkbook.Worksheets.Count))
        GetSheet.Name = nameIn
    End If
End Function

Private Sub WriteIteration(ByVal dimN As Long, ByVal it As Long, ByVal fos As Double, ByVal delta As Double, ByVal move As Double, ByVal projections As Long)
    Dim ws As Worksheet, r As Long
    Set ws = GetSheet("VBA Iterations")
    r = ws.Cells(ws.Rows.Count, 1).End(xlUp).Row + 1
    ws.Cells(r, 1).Resize(1, 6).Value = Array(dimN, it, fos, delta, move, projections)
End Sub

Private Sub WriteElements(ByVal sheetName As String, ByRef els() As ModelElement, ByRef result As ModelResult, ByRef cfg As ModelConfig)
    Dim ws As Worksheet, arr() As Variant, i As Long, n As Long
    Set ws = GetSheet(sheetName): ws.Cells.Clear
    ws.Range("A1").Value2 = "Primary VBA states; x/y/z m, area m2, forces kN, normal/c/shear kPa, principal stresses MPa."
    ws.Range("A2").Value2 = result.status & ": " & result.reason
    ws.Range("A4:W4").Value = Array("x", "y", "z_base", "height", "A_base", "W", "u", "g", "kappa", "rho", "trial normal", _
        "actual final normal", "phi deg", "c kPa", "tau kPa", "sigma1 MPa", "sigma2 MPa", "sigma3 MPa", "angle rad", "EJC MPa", "closure MPa", "eq/W", "status")
    n = result.active
    If n < 1 Then
        ws.Range("A5").Value2 = "NOT EVALUATED: no active elements for this run."
        ws.Columns("A:W").AutoFit
        Exit Sub
    End If
    ReDim arr(1 To n, 1 To 23)
    For i = 1 To n
        With els(i)
            arr(i, 1) = .x: arr(i, 2) = .y: arr(i, 3) = .zb: arr(i, 4) = .h: arr(i, 5) = .area
            arr(i, 6) = .weight: arr(i, 7) = .pressure: arr(i, 8) = .g: arr(i, 9) = .kap: arr(i, 10) = .rho
            arr(i, 11) = .sn: arr(i, 13) = .phi: arr(i, 14) = .cohesion: arr(i, 23) = .label
            If result.evaluated Then
                arr(i, 12) = .finalN: arr(i, 15) = .tau: arr(i, 22) = .equilibrium
                If cfg.strength = 1 Then
                    arr(i, 16) = .s1: arr(i, 17) = .s2: arr(i, 18) = .s3
                    arr(i, 19) = .angular: arr(i, 20) = .ejc: arr(i, 21) = .closureRes
                End If
            Else
                arr(i, 12) = "NOT EVALUATED": arr(i, 15) = "NOT EVALUATED"
            End If
        End With
    Next i
    ws.Range("A5").Resize(n, 23).Value = arr: ws.Columns("A:W").AutoFit
End Sub

Private Sub WriteSummary(ByVal ws As Worksheet, ByVal firstRow As Long, ByRef result As ModelResult)
    ws.Cells(firstRow, 2).Value2 = result.status
    ws.Cells(firstRow + 1, 4).Value2 = "WITHHELD"
    If result.status = "ACCEPTED" Then ws.Cells(firstRow + 1, 4).Value2 = result.fos
    ws.Cells(firstRow + 1, 2).Formula = "=IFERROR(IF(C8=""CURRENT"",D" & (firstRow + 1) & ",""WITHHELD""),""WITHHELD"")"
    ws.Cells(firstRow + 2, 2).Value2 = result.active
    ws.Cells(firstRow + 3, 2).Value2 = result.reason
    ws.Cells(firstRow + 4, 2).Value2 = "NOT EVALUATED"
    If result.evaluated Then ws.Cells(firstRow + 4, 2).Value2 = result.tension & " tension elements; min " & Format(result.minNormal, "0.000") & " kPa"
End Sub

Private Sub InvalidateOutputStates(ByVal explanation As String)
    Dim nameIn As Variant, ws As Worksheet
    On Error Resume Next
    For Each nameIn In Array("VBA Columns", "VBA Slices")
        Set ws = Nothing
        Set ws = ThisWorkbook.Worksheets(CStr(nameIn))
        If Not ws Is Nothing Then
            ws.Cells.Clear
            ws.Range("A1").Value2 = "NOT EVALUATED: no current primary state table."
            ws.Range("A2").Value2 = explanation
        End If
    Next nameIn
End Sub

Private Sub WriteInputRecord(ByVal mesh3 As Long, ByVal mesh2 As Long, ByVal maxIter As Long)
    Dim ws As Worksheet, inputs As Worksheet, addresses As Variant, i As Long
    Set ws = GetSheet("VBA Input Record"): Set inputs = ThisWorkbook.Worksheets("Interface")
    ws.Cells.Clear
    ws.Range("A1").Value2 = "Frozen inputs for the last primary VBA run"
    ws.Range("B1").Value2 = Format(Now, "yyyy-mm-dd hh:nn:ss")
    ws.Range("A4").Value2 = "3D mesh": ws.Range("B4").Value2 = mesh3
    ws.Range("A5").Value2 = "2D mesh": ws.Range("B5").Value2 = mesh2
    ws.Range("A6").Value2 = "Iteration limit": ws.Range("B6").Value2 = maxIter
    ws.Range("A9:C9").Value = Array("Interface address", "frozen value", "formula/literal text")
    addresses = Array("C6", "D6", "E6", "C9", "D14", "D16", "E16", "D20", "D21", "D22", "D23", "D24", "D25", _
        "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P14", "P16", "P17", "P18", "P19")
    For i = LBound(addresses) To UBound(addresses)
        ws.Cells(i + 10, 1).Value2 = addresses(i)
        ws.Cells(i + 10, 2).Value2 = inputs.Range(addresses(i)).Value2
        ws.Cells(i + 10, 3).Value2 = "'" & CStr(inputs.Range(addresses(i)).Formula)
    Next i
    ws.Columns("A:C").AutoFit
End Sub

Public Sub Thesis_RunVBAEngine()
    Dim cfg As ModelConfig, els3() As ModelElement, els2() As ModelElement, r3 As ModelResult, r2 As ModelResult
    Dim ws As Worksheet, logWs As Worksheet, n3 As Long, n2 As Long, oldScreen As Boolean, errText As String
    On Error GoTo Failed
    oldScreen = Application.ScreenUpdating: Application.ScreenUpdating = False
    Application.CalculateFullRebuild
    Set ws = GetSheet("VBA Run")
    ws.Range("B8").Value2 = "RUNNING"
    ws.Range("B11").Value2 = "WITHHELD": ws.Range("B18").Value2 = "WITHHELD"
    ws.Range("B12,B14,B19,B21").Value2 = "NOT EVALUATED"
    ws.Range("B13,B20").Value2 = "RUNNING"
    InvalidateOutputStates "A new primary run is in progress."
    ws.Range("B26").Value2 = "NOT EVALUATED": ws.Range("D26").Value2 = "NOT EVALUATED"
    ReadConfig cfg
    n3 = IntegerAt(ws, "B4"): n2 = IntegerAt(ws, "B5"): cfg.maxIter = IntegerAt(ws, "B6")
    WriteInputRecord n3, n2, cfg.maxIter
    If cfg.maxIter < 1 Or cfg.maxIter > 50000 Then Err.Raise 5, , "Iteration limit must be 1-50000."
    Set logWs = GetSheet("VBA Iterations")
    logWs.Cells.Clear
    logWs.Range("A1:F1").Value = Array("dimension", "iteration", "F", "dF/F", "max normal move kPa", "path projections")
    ws.Range("B8").Value2 = "RUNNING - Excel " & Application.Version
    SolveModel cfg, 3, n3, els3, r3, True: SolveModel cfg, 2, n2, els2, r2, True
    WriteSummary ws, 10, r3: WriteSummary ws, 17, r2
    WriteElements "VBA Columns", els3, r3, cfg: WriteElements "VBA Slices", els2, r2, cfg
    ws.Range("B8").Value2 = "RUN COMPLETE - Excel " & Application.Version & " - " & Format(Now, "yyyy-mm-dd hh:nn:ss")
    ws.Range("B25").Value2 = "No result is a thesis-closure FoS unless P4=1, P5=1, P19=1 and that dimension is ACCEPTED."
    ws.Range("D26").Value2 = "NOT EVALUATED"
    If r3.status = "ACCEPTED" Then ws.Range("D26").Value2 = r3.moment
    ws.Range("B26").Formula = "=IFERROR(IF(AND(C8=""CURRENT"",B10=""ACCEPTED""),D26,""NOT EVALUATED""),""NOT EVALUATED"")"
    ws.Range("C8").Calculate: ws.Range("B11,B18,B26").Calculate
    GoTo Cleanup
Failed:
    errText = "Error " & Err.Number & ": " & Err.Description
    On Error Resume Next
    InvalidateOutputStates "The new run failed: " & errText
    If Not ws Is Nothing Then
        ws.Range("B8").Value2 = "FAILED - " & errText
        ws.Range("B10").Value2 = "NOT ACCEPTED": ws.Range("B17").Value2 = "NOT ACCEPTED"
        ws.Range("B11").Value2 = "WITHHELD": ws.Range("B18").Value2 = "WITHHELD"
        ws.Range("D11").Value2 = "WITHHELD": ws.Range("D18").Value2 = "WITHHELD"
        ws.Range("B12,B14,B19,B21").Value2 = "NOT EVALUATED"
        ws.Range("B13,B20").Value2 = errText
        ws.Range("B26").Value2 = "NOT EVALUATED": ws.Range("D26").Value2 = "NOT EVALUATED"
    End If
Cleanup:
    Application.StatusBar = False: Application.ScreenUpdating = oldScreen
    If Len(errText) > 0 Then MsgBox errText, vbExclamation
End Sub

Private Sub BenchConfig(ByRef cfg As ModelConfig)
    cfg.xc = 355.7423115464064: cfg.zc = 644.4192948645371: cfg.radius = 725.2669975187362
    cfg.height = 360: cfg.slope = Tan(50 * PI_C / 180): cfg.gamma = 27
    cfg.sci = 77.7: cfg.mb = 10 * Exp(-58# / 28): cfg.s = Exp(-58# / 9): cfg.a = 0.5: cfg.nu = 0.271
    cfg.tolA = 0.0000000001: cfg.tolE = 0.00000001: cfg.tolF = 0.00000001
    cfg.strength = 1: cfg.closure = 3: cfg.checks = 1: cfg.maxIter = 2000
    cfg.uExtra = 0: cfg.waterAngle = 0: cfg.waterHeight = 0: cfg.cMC = 100: cfg.phiMC = 0
End Sub

Public Sub Thesis_RunEngineQA()
    Dim ws As Worksheet, cfg As ModelConfig, els3() As ModelElement, els2() As ModelElement
    Dim r3 As ModelResult, r2 As ModelResult, k As Long, passN As Long, nameIn As String
    Dim e3 As Double, e2 As Double, ok As Boolean, rowN As Long
    On Error GoTo Failed
    Set ws = GetSheet("Native Engine QA"): ws.Cells.Clear
    ws.Range("A1").Value2 = "v14 primary VBA assembly regressions - Excel " & Application.Version & " - " & Format(Now, "yyyy-mm-dd hh:nn:ss")
    ws.Range("A3:H3").Value = Array("case", "expected 3D", "actual 3D", "expected 2D", "actual 2D", "3D status", "2D status", "result")
    For k = 1 To 11
        BenchConfig cfg: e3 = 0: e2 = 0
        r3 = EmptyResult(): r2 = EmptyResult()
        Select Case k
            Case 1: nameIn = "Thesis closure Chile: no admissible FoS": cfg.closure = 1
            Case 2: nameIn = "HB dry a=.5": e3 = 2.5536684200376: e2 = 2.31025502757599
            Case 3: nameIn = "HB wet a=.5": cfg.waterAngle = 50: cfg.waterHeight = 360
                e3 = 1.20012184467095: e2 = 1.03587083245705
            Case 4: nameIn = "HB dry a=.509923": cfg.a = 0.509923: e3 = 2.47603808311863: e2 = 2.24206326532995
            Case 5, 9: nameIn = "HB wet a=.509923": cfg.a = 0.509923: cfg.waterAngle = 50: cfg.waterHeight = 360
                e3 = 1.15842745018062: e2 = 1.00046615780691
            Case 6: nameIn = "Published circle": cfg.xc = 207.28: cfg.zc = 586.53: cfg.radius = 622.08
                e3 = 2.14193739499583: e2 = 2.00387341089056
            Case 7: nameIn = "MC zero friction": cfg.strength = 0: e3 = 0.255354533635107: e2 = 0.190263799071471
            Case 8: nameIn = "Moscarda inadmissibility: no accepted FoS": cfg.closure = 2
            Case 10: nameIn = "Original MC geometry: refused, no FoS"
                cfg.xc = -1.45: cfg.zc = 3.24: cfg.radius = 4.46: cfg.height = 2.1: cfg.slope = 1 / 3.5
                cfg.gamma = 19: cfg.strength = 0: cfg.cMC = 7.2: cfg.phiMC = 30
            Case 11: nameIn = "Positive-centre cropped geometry: refused, no FoS"
                cfg.xc = 207.28: cfg.zc = 300: cfg.radius = 622.08
        End Select
        SolveModel cfg, 3, 20, els3, r3: SolveModel cfg, 2, 40, els2, r2
        If e3 > 0 Then
            ok = (r3.status = "ACCEPTED" And r2.status = "ACCEPTED")
            ok = ok And Abs(r3.fos - e3) <= 0.000000001 * e3 And Abs(r2.fos - e2) <= 0.000000001 * e2
        Else
            ok = (r3.status = "NOT ACCEPTED" And r2.status = "NOT ACCEPTED")
        End If
        If k >= 10 Then
            ok = ok And r3.active = 0 And r2.active = 0 And Not r3.evaluated And Not r2.evaluated
            ok = ok And r3.fos = 0 And r2.fos = 0 And Left$(r3.reason, 21) = "UNSUPPORTED GEOMETRY:" And Left$(r2.reason, 21) = "UNSUPPORTED GEOMETRY:"
        End If
        If k = 9 Then
            nameIn = "Final VBA normal +1 kPa: rejected"
            els3(1).finalN = els3(1).finalN + 1
            ok = ok And Not ValidateFinal(cfg, els3, r3.active, r3.fos, r3, False)
            r3.status = "NOT ACCEPTED (PERTURBED)": e3 = 0
        End If
        rowN = k + 3
        ws.Cells(rowN, 1).Resize(1, 8).Value = Array(nameIn, IIf(e3 > 0, e3, "WITHHELD"), IIf(r3.status = "ACCEPTED", r3.fos, "WITHHELD"), _
            IIf(e2 > 0, e2, "WITHHELD"), IIf(r2.status = "ACCEPTED", r2.fos, "WITHHELD"), r3.status, r2.status, IIf(ok, "PASS", "FAIL"))
        If ok Then passN = passN + 1
    Next k
    ws.Range("A15").Value2 = passN & " of 11 primary engine checks": ws.Range("H15").Value2 = IIf(passN = 11, "ALL PASS", "FAIL")
    ws.Columns("A:H").AutoFit: Application.StatusBar = False
    Exit Sub
Failed:
    On Error Resume Next
    ws.Range("H15").Value2 = "FAIL - Error " & Err.Number & ": " & Err.Description
    Application.StatusBar = False
End Sub

Private Function EmptyResult() As ModelResult
    Dim r As ModelResult
    EmptyResult = r
End Function

Public Sub Thesis_CompletePackageQA()
    Dim ws As Worksheet, formulaQA As Worksheet, engineQA As Worksheet
    On Error GoTo Failed
    Set ws = GetSheet("Completion Status")
    ws.Range("B4").Value2 = "RUNNING"
    ws.Range("B21").Value2 = "v14: Thesis_CompletePackageQA. Inspect 22 strength fixtures, 10 worksheet checks and 11 primary-engine checks."
    Thesis_RunFinalQA
    Thesis_RunEngineQA
    Set formulaQA = GetSheet("Native Workbook QA"): Set engineQA = GetSheet("Native Engine QA")
    If formulaQA.Range("L15").Value2 = "ALL PASS" And engineQA.Range("H15").Value2 = "ALL PASS" Then
        ws.Range("B4").Value2 = "NATIVE NUMERICAL CHECKS PASS"
    Else
        ws.Range("B4").Value2 = "FAIL - inspect native QA sheets"
    End If
    ws.Range("B5").Value2 = "Excel " & Application.Version & "; " & Format(Now, "yyyy-mm-dd hh:nn:ss")
    Exit Sub
Failed:
    On Error Resume Next
    ws.Range("B4").Value2 = "FAIL - " & Err.Number & ": " & Err.Description
End Sub
