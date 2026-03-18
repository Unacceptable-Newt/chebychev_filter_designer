mod cheb_calc;
mod attenuation;
mod complex;

use crate::cheb_calc::*;
use crate::attenuation::*;

use std::io;

use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyEventKind};
use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Stylize, Color},
    symbols::{border, Marker},
    layout::{Constraint, Layout, Direction},
    text::{Line, Text, Span},
    widgets::{Block, Paragraph, Widget, Dataset, GraphType, Axis, Chart},
    DefaultTerminal, Frame,
};

#[derive(Debug)]
pub struct App{
    filt_start_freq: f64,
    filt_end_freq: f64,
    filt_stages: usize,
    filt_ripple: f64,
    graph_points: usize,
    graph_start_freq: f64,
    graph_end_freq: f64,
    option_select: usize,
    incr_mult: f64,
    exit: bool,
}

impl Default for App {
    fn default() -> Self {
        Self{
            filt_start_freq: 1500000f64,
            filt_end_freq:   3500000f64,
            filt_stages: 5,
            filt_ripple: 0.5f64,
            graph_points: 1000,
            graph_start_freq: 5000f64,
            graph_end_freq: 5000000f64,
            option_select: 0,
            incr_mult: 10000.00f64,
            exit: false,
        }
    }
}

impl App {
    pub fn run(&mut self, terminal: &mut DefaultTerminal) -> io::Result<()> {
        while !self.exit {
            terminal.draw(|frame| self.draw(frame))?;
            self.handle_events()?;
        }
        Ok(())
    }

    pub fn draw(&self, frame: &mut Frame) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints(vec![
                Constraint::Fill(1),
                Constraint::Length(3),
            ]).split(frame.area());
        frame.render_widget(self, layout[1]);
        self.render_graph(frame, layout[0]);
    }

    fn handle_events(&mut self) -> io::Result<()> {
        match event::read()? {
            Event::Key(key_event) if key_event.kind == KeyEventKind::Press => {
                self.handle_key_event(key_event)
            }
            _ => {}
        };
        Ok(())
    }

    fn handle_key_event(&mut self, key_event: KeyEvent) {
        match key_event.code {
            KeyCode::Char('q') => self.exit(),
            KeyCode::Left => self.change_left(),
            KeyCode::Right => self.change_right(),
            KeyCode::Up => self.change_up(),
            KeyCode::Down => self.change_down(),
            KeyCode::Char('m') => self.incr_incr_mult(),
            KeyCode::Char('n') => self.dec_incr_mult(),
            _ => {}
        }
    }

    fn exit(&mut self) {
        self.exit = true;
    }

    fn change_left(&mut self) {
        self.option_select = 
            if self.option_select <= 0 {
                0
            }
            else {
                self.option_select - 1
            };
    }

    fn change_right(&mut self) {
        self.option_select =
            if self.option_select < 3 {
                self.option_select + 1
            }
            else {
                self.option_select
            };
    }

    fn change_up(&mut self) {
        match self.option_select {
            0 => self.inc_start_freq(),
            1 => self.inc_end_freq(),
            2 => self.inc_stages(),
            3 => self.inc_ripple(),
            _ => {}
        }
    }

    fn change_down(&mut self) {
        match self.option_select {
            0 => self.dec_start_freq(),
            1 => self.dec_end_freq(),
            2 => self.dec_stages(),
            3 => self.dec_ripple(),
            _ => {}
        }
    }

    fn inc_start_freq(&mut self) {
        self.filt_start_freq += self.incr_mult;
    }

    fn dec_start_freq(&mut self) {
        self.filt_start_freq -= self.incr_mult;
    }
    fn inc_end_freq(&mut self) {
        self.filt_end_freq += self.incr_mult;
    }

    fn dec_end_freq(&mut self) {
        self.filt_end_freq -= self.incr_mult;
    }

    fn inc_stages(&mut self) {
        self.filt_stages += 1;
    }

    fn dec_stages(&mut self) {
        self.filt_stages -= 1;
    }

    fn inc_ripple(&mut self) {
        self.filt_ripple += 0.01;
    }

    fn dec_ripple(&mut self) {
        self.filt_ripple -= 0.01;
    }

    fn incr_incr_mult(&mut self) {
        self.incr_mult *= 10f64;
    }

    fn dec_incr_mult(&mut self) {
        self.incr_mult /= 10f64;
    }

    fn render_graph(&self, frame:&mut Frame, area: Rect) {
        let (_bp_gs, bp_lcs) = chebyshev_bandpass_elements(
            self.filt_stages,
            self.filt_ripple,
            self.filt_start_freq,
            self.filt_end_freq,
            );
        let attns: Vec<(f64,f64)> = (0..self.graph_points)
            .map(|x| {
                let x = self.graph_start_freq + (self.graph_end_freq - self.graph_start_freq) * (x as f64) / (self.graph_points as f64);
                (x,-attenuation_db(&bp_lcs,x,50f64))
            }
            ).collect();
        let ds = Dataset::default()
            .name("ideal filter")
            .marker(Marker::Braille)
            .graph_type(GraphType::Line)
            .style(Color::Blue)
            .data(&attns);
        let x_axis = Axis::default()
            .title("Frequency")
            .bounds([self.graph_start_freq,self.graph_end_freq])
            .labels([format!("{:.02}",self.graph_start_freq),format!("{:.02}",self.graph_end_freq)]);
        let y_axis = Axis::default()
            .title("Attenuation")
            .bounds([-10f64,0f64])
            .labels(["-10f64","0f64"]);

        let chart = Chart::new(vec![ds]).x_axis(x_axis).y_axis(y_axis);
        frame.render_widget(chart, area);
    }
}

impl Widget for &App {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let instructions = Line::from(vec![
            " Select L ".into(),
            "<Left>".blue().bold(),
            " Select R ".into(),
            "<Right>".blue().bold(),
            " Incr Selected ".into(),
            "<Up>".blue().bold(),
            " Dec Selected ".into(),
            "<Down>".blue().bold(),
            " Increce Mag ".into(),
            "<n>".blue().bold(),
            " Decreace Mag ".into(),
            "<m>".into(),
            " Quit ".into(),
            "<Q> ".blue().bold(),
        ]);
        let block = Block::bordered()
            .title_bottom(instructions.centered())
            .border_set(border::THICK);

        let mut options: Vec<Span>= vec![
            "FS_freq: ".into(),
            format!("{:.2}",self.filt_start_freq).yellow(),
            " FE_freq: ".into(),
            format!("{:.2}",self.filt_end_freq).yellow(),
            " n_filt: ".into(),
            self.filt_stages.to_string().yellow(),
            " ripple_filt: ".into(),
            format!("{:.2}",self.filt_ripple).yellow(),
            " increace order: ".into(),
            format!("{:.2}",self.incr_mult).yellow(),
        ];

        for (i, o) in options.iter_mut().enumerate() {
            if i == (self.option_select * 2) {
                *o = o.clone().bold();
            }
        }


        let counter_text = Text::from(vec![Line::from(options)]);

        Paragraph::new(counter_text)
            .centered()
            .block(block)
            .render(area, buf);

    } 
}

#[cfg(test)]
mod tests {
    use super::*;
    use ratatui::style::Style;

    #[test]
    fn render() {
        let app = App::default();
        let mut buf = Buffer::empty(Rect::new(0, 0, 50, 4));

        app.render(buf.area, &mut buf);

        let mut expected = Buffer::with_lines(vec![
            "┏━━━━━━━━━━━━━ Counter App Tutorial ━━━━━━━━━━━━━┓",
            "┃                    Value: 0                    ┃",
            "┃                                                ┃",
            "┗━ Decrement <Left> Increment <Right> Quit <Q> ━━┛",
        ]);
        let title_style = Style::new().bold();
        let counter_style = Style::new().yellow();
        let key_style = Style::new().blue().bold();
        expected.set_style(Rect::new(14, 0, 22, 1), title_style);
        expected.set_style(Rect::new(28, 1, 1, 1), counter_style);
        expected.set_style(Rect::new(13, 3, 6, 1), key_style);
        expected.set_style(Rect::new(30, 3, 7, 1), key_style);
        expected.set_style(Rect::new(43, 3, 4, 1), key_style);

        assert_eq!(buf, expected);
    }

    #[test]
    fn handle_key_event() {
        let mut app = App::default();
        app.handle_key_event(KeyCode::Right.into());
        assert_eq!(app.counter, 1);

        app.handle_key_event(KeyCode::Left.into());
        assert_eq!(app.counter, 0);

        let mut app = App::default();
        app.handle_key_event(KeyCode::Char('q').into());
        assert!(app.exit);
    }
}

fn main() -> io::Result<()>{
    ratatui::run(|terminal| App::default().run(terminal))
}
