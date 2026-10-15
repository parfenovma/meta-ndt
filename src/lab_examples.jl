export lab_examples, example_experiment, example_reference

"""Six frozen specimens; questions describe experiments, not guaranteed metamaterial effects."""
function lab_examples()
    [
        (id=:single_cylinder,title="1 · Один цилиндр",question="Что меняет упругость?",
         description="Пустая вода, жёсткое препятствие и движущийся алюминий — базовый контроль сопряжения."),
        (id=:cylinder_lattice,title="2 · Решётка цилиндров",question="Как элементы взаимодействуют?",
         description="Девять включений решаются совместно. Смотрите поле между ними и частотную зависимость сигнала."),
        (id=:missing_cell,title="3 · Пропущенная ячейка",question="Как меняется поле при дефекте решётки?",
         description="Из той же решётки удалён центральный цилиндр. Дополнительный контроль — полная решётка; локализация не предполагается заранее."),
        (id=:straight_brick,title="4 · Прямой кирпичик",question="Как выглядит контроль для STG?",
         description="Длина и площадь совпадают с одиночным STG: сравнение формы при одинаковой массе материала."),
        (id=:stg_single,title="5 · Синусоидальный STG",question="Как извилистость меняет колебания и рассеяние?",
         description="Настоящие сдвинутые профили: L=17, H=7, A=2.5 mm, N=2; H−A=4.5 mm в общей части. Контроль — прямой кирпичик равной массы."),
        (id=:stg_row,title="6 · Ряд STG-кирпичиков",question="Что меняет взаимодействие трёх кирпичиков?",
         description="Три отдельных упругих STG в воде. Контроль — три прямых кирпичика той же массы, длины и положения."),
    ]
end

function example_experiment(id::Symbol; frequency=100e3)
    frequency >= 40e3 && frequency <= 160e3 || throw(ArgumentError("example carrier must be within 40–160 kHz"))
    centers=[(0.0,0.0)]; shape=nothing
    if id in (:cylinder_lattice,:missing_cell)
        centers=cylinder_lattice(3,3,.018)
        id == :missing_cell && filter!(p -> p != (0.0,0.0),centers)
    elseif id in (:straight_brick,:stg_single,:stg_row)
        shape=SinusoidalBrick()
        if id == :straight_brick
            shape=RectangularBrick(length=shape.length,height=specimen_area(shape)/shape.length)
        elseif id == :stg_row
            centers=[(0.0,-.018),(0.0,0.0),(0.0,.018)]
        end
    elseif id != :single_cylinder
        throw(ArgumentError("unknown laboratory example"))
    end
    s=CylinderExperiment(radius=.004,centers=centers,shape=shape,frequency=frequency,
        source=(-.055,0.0),source_width=.002,
        sensors=[(.040,0.0),(-.040,.020),(.045,-.025)],
        mesh_size=min(.001,1480/(14frequency)),time_step=min(.25e-6,.04/frequency))
    validate_experiment(s)
end

function example_reference(id::Symbol,s::CylinderExperiment)
    id in (:stg_single,:stg_row) && return straight_reference(s)
    id == :missing_cell && return remake_experiment(s;centers=cylinder_lattice(3,3,.018))
    nothing
end
